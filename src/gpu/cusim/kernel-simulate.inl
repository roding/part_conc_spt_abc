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
		class tSimulateSetup,
		class tSimulateRun,
		class tOutput,
		class tRandData
	> 
	__global__
	void /*__launch_bounds__(256,8)*/ K_simulate_system( tSimulateSetup aSetup, tSimulateRun aRunData, tOutput aOutput, tRandData aRandData )
	{
		using Real_ = typename tSimulateSetup::Scalar;
		using Count_ = typename tSimulateSetup::Count;

		static_assert( std::is_same<Real_,typename tRandData::result_type>::value, "tSimulateSetup::Scalar and random number generator result type mismatch" );
		static_assert( std::is_same<Count_,typename tOutput::value_type>::value, "tSimulateSetup::Count and output type mismatch" );

		constexpr unsigned kWarpSize = 32u; /* WARNING: warp size assumption */
		
		unsigned const gtid = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned const warp = gtid / kWarpSize; /* WARNING: layout assumption */
		unsigned const wtid = gtid % kWarpSize; /* WARNING: layout assumption */

		if( warp >= aSetup.jobCount )
			return;

		using RandomFactory_ = typename tRandData::Factory;
		auto ran = typename RandomFactory_::Instance( gtid, aRandData );

		auto const frames = aRunData.frames[warp];
		auto const particles = aRunData.particles[warp];

		for( Count_ particle = wtid; particle < particles; particle += kWarpSize )
		{
			auto part = make_particle( aSetup, aRunData, gtid, ran, aRandData );
			
			Real_ de = Real_(0);
			Count_ detectionCount = detected( part, aSetup, aRunData )
				? Count_(1)
				: Count_(0)
			;

			for( Count_ frame = 1; frame < frames; ++frame )
			{
				auto const rwsd = rand_walk_stddev( part, aSetup, aRunData );

				auto const dx = rwsd * ran.normal01( gtid, aRandData );
				part.x = periodic_bounds_( part.x + dx, aSetup.halfLx );
				
				auto const dy = rwsd * ran.normal01( gtid, aRandData );
				part.y = periodic_bounds_( part.y + dy, aSetup.halfLy );

				auto const dz = rwsd * ran.normal01( gtid, aRandData );
				part.z = periodic_bounds_( part.z + dz, aSetup.halfLz );

				if( detected( part, aSetup, aRunData ) )
				{
					++detectionCount;

					if( detectionCount >= Count_(2) )
						de += dx*dx + dy*dy;
				}
				else
				{
					if( detectionCount >= aSetup.kmin )
					{
						aOutput.record( detectionCount, de / ((detectionCount-1) * Real_(4) * aSetup.deltaT) );
					}

					detectionCount = Count_(0);
					de = Real_(0);
				}
			}

			if( detectionCount >= aSetup.kmin )
			{
				aOutput.record( detectionCount, de / ((detectionCount-1) * Real_(4) * aSetup.deltaT) );
			}
		}

		ran.store( gtid, aRandData );
	}
}

