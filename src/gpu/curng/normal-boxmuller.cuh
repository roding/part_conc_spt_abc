/*-******************************************************* -- HEADER -{{{1- */
/*-	Normal distribution: Box-Muller
 */
/*-***************************************************************** -}}}1- */

#ifndef NORMAL_BOXMULLER_CUH_BB2157D5_F832_4771_BA4F_C71C0E067A69
#define NORMAL_BOXMULLER_CUH_BB2157D5_F832_4771_BA4F_C71C0E067A69

#include "support.cuh"
#include "../support/utility.hpp"

namespace curng
{
	/** Normal distribution: Box-Muller
	 * 
	 * Generates normally distributed numbers with mean 0 and an unit standard 
	 * deviation using the Box-Muller method.
	 *
	 * Usage:
	   \code
	   // host
       auto normalData = curng::initialize( 
	       Identity<NormalBoxMuller>,
		   Identity<float>,
		   numberOfCUDAThreads
	   );

	   ...
	   kernel<<<...>>>( normalData, ... );
	   ...

	   cleanup( normalData );

	   // device
	   template< class tDist, class tEngine >
	   __global__ void kernel( tDist& aDistData, tEngine& aEngineData )
	   {
		   auto tid = threadIdx.x + blockIdx.x*blockDim.x + ...; // global thread ID
           auto engine = typename tEngine::Engine( tid, aEngineData );
           auto normalDist = typename tDist::Distribution<float>( tid, aDistData );

		   ...
           float number = normalDist( engine, tid, aDistData, aEngineData );
		   ...
	   }
	   \endcode
	 *
	 * \note The methods uses a single `tReal` element cache; two numbers are
	 * generated on every second call to `operator()`. One is returned
	 * immediately, and the other one is returned on the next call to
	 * `operator()`. This is a consequence of the Box-Muller method, which
	 * always yields two independent normally distributed numbers.
	 */
	struct NormalBoxMuller
	{
		template< typename tReal >
		struct GlobalData
		{
			using Factory = NormalBoxMuller;
			using result_type = tReal;
		};

		template< typename tReal >
		class Distribution
		{
			public:
				using result_type = tReal;

			public:
				__device__
				Distribution( unsigned, GlobalData<tReal> const& );

			public:
				template< class tRand, class tRandData > __device__
				result_type operator() (tRand&, unsigned, GlobalData<tReal>&, tRandData&);

			private:
				tReal mCache;
		};

		template< typename tReal > static __host__
		GlobalData<tReal> initialize( Identity<tReal>, std::size_t );

		template< typename tReal > static __host__
		void cleanup( GlobalData<tReal>& );
	};
}

#include "normal-boxmuller.inl"
#endif // NORMAL_BOXMULLER_CUH_BB2157D5_F832_4771_BA4F_C71C0E067A69


