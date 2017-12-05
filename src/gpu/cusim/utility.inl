#include "../support/dsvalue.hpp"

namespace cusim
{
	// This could would be improved vastly by "if constexpr". But, alas, the
	// future is not here, yet.
	
	// wrand_index()
	namespace detail
	{
		template< typename tReal, class tCompCount >
		struct WrandIndexImpl_
		{
			static __device__
			unsigned enact( tReal aValue, tCompCount aCount, DSInline<tReal,tCompCount> const& aPreWeights )
			{
				unsigned bot = 0;
				unsigned top = aCount;

				// Binary search
				while( bot != top )
				{
					unsigned const mid = (top+bot)/2;
					if( aPreWeights.values[mid] <= aValue )
						bot = mid+1;
					else
						top = mid;
				}
				
				return bot;
			}
		};

		template< typename tReal >
		struct WrandIndexImpl_<tReal, StaticValue<unsigned,2>>
		{
			static __device__
			unsigned enact( tReal aValue, StaticValue<unsigned,2>, DSInline<tReal,StaticValue<unsigned,2>> const& aPreWeights )
			{
				if( aPreWeights.values[0] < aValue )
					return 0;
				
				return 1;
			}
		};

		// We should never call this with tCompCount = 1, since doing so would
		// likely imply that we generated a random number that's not neccessary.
		template< typename tReal >
		struct WrandIndexImpl_<tReal, StaticValue<unsigned,1>>;
	}

	template< typename tReal, class tCompCount > __device__ inline
	unsigned wrand_index( tReal aValue, tCompCount aCount, DSInline<tReal,tCompCount> const& aPreWeights )
	{
		return detail::WrandIndexImpl_<tReal,tCompCount>::enact( aValue, aCount, aPreWeights );
	}


	// make_particle()
	namespace detail
	{
		template< class tCompCount >
		struct MakeParticleImpl_
		{
			template< class tSystem, class tSimRun, class tRand, class tRandData > static __device__
			auto enact( tSystem const& aSystem, tSimRun const& aRunData, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSystem::Particle
			{
				typename tSystem::Particle ret;
				ret.x = aSystem.halfLx * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.y = aSystem.halfLy * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.z = aSystem.halfLz * (2*aRan.uniform01( aTid, aRanData )-1);

				auto const rr = aRan.uniform01( aTid, aRanData );
				ret.index = wrand_index( rr, aSystem.compCount, aRunData.preCompProb );
				return ret;
			}
		};

		template<>
		struct MakeParticleImpl_< StaticValue<unsigned,1> >
		{
			template< class tSystem, class tSimRun, class tRand, class tRandData > static __device__
			auto enact( tSystem const& aSystem, tSimRun const&, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSystem::Particle
			{
				typename tSystem::Particle ret;
				ret.x = aSystem.halfLx * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.y = aSystem.halfLy * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.z = aSystem.halfLz * (2*aRan.uniform01( aTid, aRanData )-1);
				return ret;
			}
		};
	}

	template< class tSimSetup, class tSimRun, class tRand, class tRandData > __device__
	auto make_particle( tSimSetup const& aSystem, tSimRun const& aRunData, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSimSetup::Particle
	{
		return detail::MakeParticleImpl_<typename tSimSetup::CompCount>::enact( aSystem, aRunData, aTid, aRan, aRanData );
	}


	// detected()
	namespace detail
	{
		template< EModel tModel, class tCompCount >
		struct DetectedImpl_;

		template< class tCompCount >
		struct DetectedImpl_< EModel::discreteFixedZ, tCompCount >
		{
			template< class tSimSetup, class tSimRun > static __device__
			bool enact( typename tSimSetup::Particle const& aParticle, tSimSetup const& aSystem, tSimRun const& aRunData )
			{
				return std::abs(aParticle.z) <= aRunData.halfAz.value
					&& std::abs(aParticle.x) <= aSystem.halfAx
					&& std::abs(aParticle.y) <= aSystem.halfAy
				;
			}
		};
		template< class tCompCount >
		struct DetectedImpl_< EModel::discreteVariableZ, tCompCount >
		{
			template< class tSimSetup, class tSimRun > static __device__
			bool enact( typename tSimSetup::Particle const& aParticle, tSimSetup const& aSystem, tSimRun const& aRunData )
			{
				return std::abs(aParticle.z) <= aRunData.halfAz.values[aParticle.index]
					&& std::abs(aParticle.x) <= aSystem.halfAx
					&& std::abs(aParticle.y) <= aSystem.halfAy
				;
			}
		};
		template<>
		struct DetectedImpl_< EModel::discreteVariableZ, StaticValue<unsigned,1> >
		{
			template< class tSimSetup, class tSimRun > static __device__
			bool enact( typename tSimSetup::Particle const& aParticle, tSimSetup const& aSystem, tSimRun const& aRunData)
			{
				return std::abs(aParticle.z) <= aRunData.halfAz.value
					&& std::abs(aParticle.x) <= aSystem.halfAx
					&& std::abs(aParticle.y) <= aSystem.halfAy
				;
			}
		};
	}

	template< class tSimSetup, class tSimRun > __device__ inline
	bool detected( typename tSimSetup::Particle const& aParticle, tSimSetup const& aSystem, tSimRun const& aRunData )
	{
		return detail::DetectedImpl_<
			tSimSetup::kModel,
			typename tSimSetup::CompCount
		>::enact( aParticle, aSystem, aRunData );
	}

	// rand_walk_stddev()
	namespace detail
	{
		template< class tCompCount >
		struct RandWalkStddevImpl_
		{
			template< class tSimSetup, class tSimRun, class tPart > static __device__
			auto enact( tPart const& aParticle, tSimSetup const&, tSimRun const& aRunData ) -> typename tSimSetup::Scalar
			{
				return aRunData.randWalkStddev.values[aParticle.index];
			}
		};

		template<>
		struct RandWalkStddevImpl_<StaticValue<unsigned,1>>
		{
			template< class tSimSetup, class tSimRun, class tPart > static __device__
			auto enact( tPart const&, tSimSetup const&, tSimRun const& aRunData ) -> typename tSimSetup::Scalar
			{
				return aRunData.randWalkStddev.value;
			}
		};
	}

	template< class tSimSetup, class tSimRun, class tPart > __device__ inline
	auto rand_walk_stddev( tPart const& aParticle, tSimSetup const& aSystem, tSimRun const& aRunData ) -> typename tSimSetup::Scalar
	{
		return detail::RandWalkStddevImpl_<
			typename tSimSetup::CompCount
		>::enact( aParticle, aSystem, aRunData );
	}
}
