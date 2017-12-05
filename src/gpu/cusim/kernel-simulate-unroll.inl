#include "utility.cuh"

/* Unrolled version of the kernel in kernel-simulate.inl.
 *
 * The frame loop is unrolled once, so that we can generate 6 normal
 * distributed random numbers via three rounds of the box-muller generation
 * up-front. The 6 numbers are consumed by the two iterations -- each iteration
 * uses three to displace x, y and z, respectively.
 * 
 * This avoids some amount of logic related to caching every other normal
 * distributed number.
 *
 * A final iteration is potentially performed at the end to ensure that we 
 * perform the correct number of iterations.
 *
 * Note that this hard-codes the normal distribution. See original code in 
 * curng/normal-boxmuller.{cuh,inl}
 */

namespace cusim
{
	namespace detail_unr
	{
		template< typename tReal > __device__ inline
		tReal periodic_bounds_( tReal aVal, tReal aMax )
		{
			// Note: assumes that a particle moves much less than 2*aMax per iteration!
			if( std::abs(aVal) > aMax )
				aVal -= tReal(2)*std::copysign( aMax, aVal );

			return aVal;
		}

		__device__ inline
		float log_( float aX )
		{
			//return logf( aX );
			return __logf( aX );
		}
		__device__ inline
		double log_( double aX )
		{
			return log( aX );
		}
		
		__device__ inline
		void sincos_(float aX, float* aSin, float* aCos)
		{
			__sincosf(aX, aSin, aCos);
		}
		__device__ inline
		void sincos_(double aX, double* aSin, double* aCos)
		{
			sincos(aX, aSin, aCos);
		}

		__device__ inline
		float sin_( float aX )
		{
			return __sinf( aX );
		}
		__device__ inline
		double sin_( double aX )
		{
			return sin( aX );
		}
	}

	template< 
		class tSimulateSetup,
		class tSimulateRun,
		class tOutput,
		class tRandData
	> 
	__global__
	void /*__launch_bounds__(256,8)*/ K_simulate_system_unr( tSimulateSetup aSetup, tSimulateRun aRunData, tOutput aOutput, tRandData aRandData )
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

			auto const rwsd = rand_walk_stddev( part, aSetup, aRunData );
		
			constexpr Real_ pi2 = Real_(2)*Real_(3.14159265357989);
			constexpr Real_ eps = std::numeric_limits<Real_>::min();

			Count_ frame = 1;
			for( ; frame < frames; frame += 2 )
			{
				Real_ n0, n1, n2, n3, n4, n5;

				{
					Real_ u1;
					do
					{
						u1 = ran.uniform01( gtid, aRandData );
					} while (u1 <= eps);

					Real_ const lu2 = std::sqrt(Real_(-2) * detail_unr::log_(u1));
					Real_ const u2 = ran.uniform01( gtid, aRandData );

					Real_ s, c;
					detail_unr::sincos_(pi2 * u2, &s, &c);

					n0 = lu2 * s;
					n1 = lu2 * c;
				}
				{
					Real_ u1;
					do
					{
						u1 = ran.uniform01( gtid, aRandData );
					} while (u1 <= eps);

					Real_ const lu2 = std::sqrt(Real_(-2) * detail_unr::log_(u1));
					Real_ const u2 = ran.uniform01( gtid, aRandData );

					Real_ s, c;
					detail_unr::sincos_(pi2 * u2, &s, &c);

					n2 = lu2 * s;
					n3 = lu2 * c;
				}
				{
					Real_ u1;
					do
					{
						u1 = ran.uniform01( gtid, aRandData );
					} while (u1 <= eps);

					Real_ const lu2 = std::sqrt(Real_(-2) * detail_unr::log_(u1));
					Real_ const u2 = ran.uniform01( gtid, aRandData );

					Real_ s, c;
					detail_unr::sincos_(pi2 * u2, &s, &c);

					n4 = lu2 * s;
					n5 = lu2 * c;
				}

				{
					auto const dx = rwsd * n0;
					part.x = detail_unr::periodic_bounds_( part.x + dx, aSetup.halfLx );
					
					auto const dy = rwsd * n1;
					part.y = detail_unr::periodic_bounds_( part.y + dy, aSetup.halfLy );

					auto const dz = rwsd * n2;
					part.z = detail_unr::periodic_bounds_( part.z + dz, aSetup.halfLz );

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
							aOutput.record( detectionCount, de / ((detectionCount-1) * aSetup.fourDeltaT) );
						}

						detectionCount = Count_(0);
						de = Real_(0);
					}
				}

				{
					auto const dx = rwsd * n3;
					part.x = detail_unr::periodic_bounds_( part.x + dx, aSetup.halfLx );
					
					auto const dy = rwsd * n4;
					part.y = detail_unr::periodic_bounds_( part.y + dy, aSetup.halfLy );

					auto const dz = rwsd * n5;
					part.z = detail_unr::periodic_bounds_( part.z + dz, aSetup.halfLz );

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
							aOutput.record( detectionCount, de / ((detectionCount-1) * aSetup.fourDeltaT) );
						}

						detectionCount = Count_(0);
						de = Real_(0);
					}
				}
			}

			if( frame < frames )
			{
				Real_ n0, n1, n2;

				{
					Real_ u1;
					do
					{
						u1 = ran.uniform01( gtid, aRandData );
					} while (u1 <= eps);

					Real_ const lu2 = std::sqrt(Real_(-2) * detail_unr::log_(u1));
					Real_ const u2 = ran.uniform01( gtid, aRandData );

					Real_ s, c;
					detail_unr::sincos_(pi2 * u2, &s, &c);

					n0 = lu2 * s;
					n1 = lu2 * c;
				}
				{
					Real_ u1;
					do
					{
						u1 = ran.uniform01( gtid, aRandData );
					} while (u1 <= eps);

					Real_ const lu2 = std::sqrt(Real_(-2) * detail_unr::log_(u1));
					Real_ const u2 = ran.uniform01( gtid, aRandData );

					Real_ s = detail_unr::sin_(pi2 * u2);
					n2 = lu2 * s;
				}

				{
					auto const dx = rwsd * n0;
					part.x = detail_unr::periodic_bounds_( part.x + dx, aSetup.halfLx );
					
					auto const dy = rwsd * n1;
					part.y = detail_unr::periodic_bounds_( part.y + dy, aSetup.halfLy );

					auto const dz = rwsd * n2;
					part.z = detail_unr::periodic_bounds_( part.z + dz, aSetup.halfLz );

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
							aOutput.record( detectionCount, de / ((detectionCount-1) * aSetup.fourDeltaT) );
						}

						detectionCount = Count_(0);
						de = Real_(0);
					}
				}
			}

			if( detectionCount >= aSetup.kmin )
			{
				aOutput.record( detectionCount, de / ((detectionCount-1) * aSetup.fourDeltaT) );
			}
		}

		ran.store( gtid, aRandData );
	}
}

