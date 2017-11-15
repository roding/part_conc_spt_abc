
namespace cusim
{
	namespace detail
	{
#		if 1
		template< typename tCount > __device__ inline
		tCount warp_prefix( unsigned aWTID, tCount aValue )
		{
#			if 0
			auto x0 = __shfl_up( aValue, 1 );
			aValue += (aWTID >= 1) ? x0 : 0;

			auto x1 = __shfl_up( aValue, 2 );
			aValue += (aWTID >= 2) ? x1 : 0;

			auto x2 = __shfl_up( aValue, 4 );
			aValue += (aWTID >= 4) ? x2 : 0;

			auto x3 = __shfl_up( aValue, 8 );
			aValue += (aWTID >= 8) ? x3 : 0;

			auto x4 = __shfl_up( aValue, 16 );
			aValue += (aWTID >= 16) ? x4 : 0;
#			else
			/* This version results in much tighter code than the C version
			 * above. The C version produces four instructions per step (a
			 * ISETP, IADD, SHFL and SEL); the inline version gets it down to
			 * two (SHFL + predicated IADD).
			 *
			 * Incidentally, the below PTX is shown in as an example for the
			 * shfl instruction in the PTX ISA document. See
			 * http://docs.nvidia.com/cuda/parallel-thread-execution/index.html
			 */
			__asm__ volatile( "{\n\t"
				".reg .u32 t0;\n\t"
				".reg .pred valid;\n\t"

				"shfl.up.b32 t0|valid, %0, 1, 0;\n\t"
				"@valid add.s32 %0, t0, %0;\n\t"

				"shfl.up.b32 t0|valid, %0, 2, 0;\n\t"
				"@valid add.s32 %0, t0, %0;\n\t"

				"shfl.up.b32 t0|valid, %0, 4, 0;\n\t"
				"@valid add.s32 %0, t0, %0;\n\t"

				"shfl.up.b32 t0|valid, %0, 8, 0;\n\t"
				"@valid add.s32 %0, t0, %0;\n\t"

				"shfl.up.b32 t0|valid, %0, 16, 0;\n\t"
				"@valid add.s32 %0, t0, %0;\n\t"

				"}\n\t"
				: "+r"(aValue)
			);
#			endif

			return aValue;
		}
#		else
		/* Early test: use shared memory to communicate values between threads.
		 * This version uses a total of 48 values per warp, with the first 16
		 * set to zero. This avoids the branches/predicated instructions seen
		 * in the code above.
		 *
		 * Brief ad-hoc benchmarking shows that it's not really worth it in
		 * terms of performance, though, as both versions clock in at roughly
		 * the same elapsed times on my test GPU.
		 *
		 * Uses inline PTX, since we want to avoid an extra fenches/syncs; this
		 * however gives the CUDA/C++ compiler freedom to mess up the method
		 * (probably rightly so).
		 */
		template< typename tCount > __device__ inline
		tCount warp_prefix( unsigned aWTID, tCount aValue, tCount* aSmem )
		{
			__asm__ volatile( "{\n\t"

				".reg .u32 t0;\n\t"
				".reg .u32 t1;\n\t"
				".reg .u32 t2;\n\t"
				".reg .u64 q0;\n\t"
				".reg .u64 q1;\n\t"
				".reg .u64 base;\n\t"

				"mov.u32 t0, %%tid.x;\n\t"
				"mul.wide.u32 q0, t0, 4;\n\t"
				"add.s64 q1, %1, q0;\n\t"
				"cvta.to.shared.u64 base, q1;\n\t"

				"st.shared.u32 [base+64], %0;\n\t"

				"ld.shared.u32 t1, [base+60];\n\t" // 64-1*4 = 60
				"add.s32 t2, %0, t1;\n\t"
				"st.shared.u32 [base+64], t2;\n\t"

				"ld.shared.u32 t1, [base+56];\n\t" // 64-2*4 = 56
				"add.s32 t2, t2, t1;\n\t"
				"st.shared.u32 [base+64], t2;\n\t"

				"ld.shared.u32 t1, [base+48];\n\t" // 64-4*4 = 48
				"add.s32 t2, t2, t1;\n\t"
				"st.shared.u32 [base+64], t2;\n\t"

				"ld.shared.u32 t1, [base+32];\n\t" // 64-8*4 = 32
				"add.s32 t2, t2, t1;\n\t"
				"st.shared.u32 [base+64], t2;\n\t"

				"ld.shared.u32 t1, [base];\n\t" // 64-16*4 = 0
				"add.s32 %0, t2, t1;\n\t"

				"}\n\t"
				: "+r"(aValue)
				: "l"(aSmem)
			);
			return aValue;
		}
#		endif

		template< typename tReal > __device__ inline
		tReal reduce0( tReal aValue )
		{
			aValue += __shfl_down( aValue, 16 );
			aValue += __shfl_down( aValue, 8 );
			aValue += __shfl_down( aValue, 4 );
			aValue += __shfl_down( aValue, 2 );
			aValue += __shfl_down( aValue, 1 );
			return aValue;
		}
	}

	template< typename tReal, typename tCount > __global__
	void /*__launch_bounds__(1024,1)*/ K_distance( tReal* aDistance, unsigned aN, unsigned aM, tCount* aCur, tCount const* aRef )
	{
		__shared__ tReal totals[32]; //XXX-FIXME: number of warps in block.
		//__shared__ tCount buff[32][48];

		auto const warp = threadIdx.y;
		auto const wtid = threadIdx.x;

		auto const n32 = (aN+32-1)/32*32;
		auto const m32 = (aM+32-1)/32*32;

		// init
		if( 0 == warp )
		{
			totals[wtid] = tReal(0);
		}

		//buff[warp][wtid] = 0;

		__syncthreads();

		// column-wise prefix sums
		for( auto row = warp; row < aN; row += blockDim.y )
		{
			tCount base = 0;
			for( auto col = wtid; col < m32; col += 32 )
			{
				tCount const val = col < aM
					? aCur[row*aM+col]
					: 0
				;

				tCount const sum = base + detail::warp_prefix( wtid, val );

				if( col < aM )
					aCur[row*aM+col] = sum;

				base = __shfl( sum, 31 );
			}
		}

		__syncthreads();

		// row-wise prefix sums, and accumulate the squared difference to the reference
		tReal acc = tReal(0);
		for( auto col = warp; col < aM; col += blockDim.y )
		{
			tCount base = 0;
			tReal a2 = tReal(0);
			for( auto row = wtid; row < n32; row += 32 )
			{
				tCount const val = row < aN
					? aCur[row*aM+col]
					: 0
				;

				tCount const sum = base + detail::warp_prefix( wtid, val );
				base = __shfl( sum, 31 );

				if( row < aN )
				{
					tCount const ref = aRef[row*aM+col];
					tReal const dd = tReal(sum) - ref;
					a2 += dd*dd;
				}
			}

			acc += a2;
		}

		// reduce the per-thread sums in each warp
		tReal const wsum = detail::reduce0( acc );
		if( 0 == wtid )
			totals[warp] = wsum;

		__syncthreads();

		// have one warp reduce the per-warp sums to the final sum
		if( 0 == warp )
		{
			tReal const tsum = wtid < 32 //XXX
				? totals[wtid]
				: 0
			;

			tReal const fin = detail::reduce0( tsum );

			if( 0 == wtid )
				*aDistance = fin;
		}

		// zero out the histogram for the next invocation
		for( auto row = warp; row < aN; row += blockDim.y )
		{
			for( auto col = wtid; col < m32; col += 32 )
			{
				aCur[row*aM+col] = 0;
			}
		}
	}
}
