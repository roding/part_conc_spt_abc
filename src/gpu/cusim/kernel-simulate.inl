#include "utility.cuh"

namespace cusim
{
	template< typename tReal > __device__ inline
	tReal periodic_bounds_( tReal aVal, tReal aMax )
	{
		// Note: assumes that a particle moves much less than 2*aMax per iteration!
		if( std::abs(aVal) > aMax ) 
			aVal -= tReal(2)*std::copysign( aMax, aVal );

		return aVal;
	}
	
	template< 
		class tSystemSetup,
		class tRunData,
		class tOutput,
		class tRandData, 
		class tReal,
		class tCount,
		class tRandomFactory
	> 
	__global__
	void /*__launch_bounds__(256,8)*/ K_simulate_system( tSystemSetup aSystem, tRunData aRunData, tOutput aOutput, tRandData aRandData )
	{
		static_assert( std::is_same<tReal,typename tRandData::result_type>::value, "System real number type and random generator real number type mismatch" );
		static_assert( std::is_same<tCount,typename tOutput::value_type>::value, "System count_type and output type mismatch" );

		constexpr unsigned kWarpSize = 32u; /* WARNING: warp size assumption */
		
		unsigned const gtid = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned const warp = gtid / kWarpSize; /* WARNING: layout assumption */
		unsigned const wtid = gtid % kWarpSize; /* WARNING: layout assumption */

		if( warp >= aSystem.jobCount )
			return;

		auto ran = typename tRandomFactory::Instance( gtid, aRandData );

		auto const frames = aRunData.frames[warp];
		auto const particles = aRunData.particles[warp];

		//TODO: maybe cache aSystem.{Lx,Ly,Lz} locally; other aSystem stuff?
		//TODO-check: register usage with/without caching

		for( tCount particle = wtid; particle < particles; particle += kWarpSize )
		{
			auto part = make_particle( aSystem, aRunData, gtid, ran, aRandData );
			
			tReal de = tReal(0);
			tCount detectionCount = detected( part, aSystem, aRunData )
				? tCount(1) 
				: tCount(0)
			;

			for( tCount frame = 1; frame < frames; ++frame )
			{
				auto const rwsd = rand_walk_stddev( part, aSystem, aRunData );

				auto const dx = rwsd * ran.normal01( gtid, aRandData );
				part.x = periodic_bounds_( part.x + dx, aSystem.halfLx );
				
				auto const dy = rwsd * ran.normal01( gtid, aRandData );
				part.y = periodic_bounds_( part.y + dy, aSystem.halfLy );

				auto const dz = rwsd * ran.normal01( gtid, aRandData );
				part.z = periodic_bounds_( part.z + dz, aSystem.halfLz );

				if( detected( part, aSystem, aRunData ) )
				{
					++detectionCount;

					if( detectionCount >= tCount(2) )
						de += dx*dx + dy*dy;
				}
				else
				{
					if( detectionCount >= aSystem.kmin )
					{
						aOutput.record( detectionCount, de / ((detectionCount-1) * tReal(4) * aSystem.deltaT) );
					}

					detectionCount = tCount(0);
					de = tReal(0);
				}
			}

			if( detectionCount >= aSystem.kmin )
			{
				aOutput.record( detectionCount, de / ((detectionCount-1) * tReal(4) * aSystem.deltaT) );
			}
		}

		ran.store( gtid, aRandData );
	}
}

