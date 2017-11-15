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
			template< class tSystem, class tRunData, class tRand, class tRandData > static __device__
			auto enact( tSystem const& aSystem, tRunData const& aRunData, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSystem::Particle
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
			template< class tSystem, class tRunData, class tRand, class tRandData > static __device__
			auto enact( tSystem const& aSystem, tRunData const&, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSystem::Particle
			{
				typename tSystem::Particle ret;
				ret.x = aSystem.halfLx * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.y = aSystem.halfLy * (2*aRan.uniform01( aTid, aRanData )-1);
				ret.z = aSystem.halfLz * (2*aRan.uniform01( aTid, aRanData )-1);
				return ret;
			}
		};
	}

	template< class tSystemSetup, class tRunData, class tRand, class tRandData > __device__
	auto make_particle( tSystemSetup const& aSystem, tRunData const& aRunData, unsigned aTid, tRand& aRan, tRandData& aRanData ) -> typename tSystemSetup::Particle
	{
		return detail::MakeParticleImpl_<typename tSystemSetup::CompCount>::enact( aSystem, aRunData, aTid, aRan, aRanData );
	}


	// detected()
	namespace detail
	{
		template< EModel tModel, class tCompCount >
		struct DetectedImpl_;

		template< class tCompCount >
		struct DetectedImpl_< EModel::discreteFixedZ, tCompCount >
		{
			template< class tSystemSetup, class tRunData > static __device__
			bool enact( typename tSystemSetup::Particle const& aParticle, tSystemSetup const& aSystem, tRunData const& aRunData )
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
			template< class tSystemSetup, class tRunData > static __device__
			bool enact( typename tSystemSetup::Particle const& aParticle, tSystemSetup const& aSystem, tRunData const& aRunData )
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
			template< class tSystemSetup, class tRunData > static __device__
			bool enact( typename tSystemSetup::Particle const& aParticle, tSystemSetup const& aSystem, tRunData const& aRunData)
			{
				return std::abs(aParticle.z) <= aRunData.halfAz.value
					&& std::abs(aParticle.x) <= aSystem.halfAx
					&& std::abs(aParticle.y) <= aSystem.halfAy
				;
			}
		};
	}

	template< class tSystemSetup, class tRunData > __device__ inline
	bool detected( typename tSystemSetup::Particle const& aParticle, tSystemSetup const& aSystem, tRunData const& aRunData )
	{
		return detail::DetectedImpl_<
			tSystemSetup::kModel,
			typename tSystemSetup::CompCount
		>::enact( aParticle, aSystem, aRunData );
	}

	// rand_walk_stddev()
	namespace detail
	{
		template< class tCompCount >
		struct RandWalkStddevImpl_
		{
			template< class tSystemSetup, class tRunData > static __device__
			auto enact( typename tSystemSetup::Particle const& aParticle, tSystemSetup const&, tRunData const& aRunData ) -> typename tSystemSetup::value_type
			{
				return aRunData.randWalkStddev.values[aParticle.index];
			}
		};

		template<>
		struct RandWalkStddevImpl_<StaticValue<unsigned,1>>
		{
			template< class tSystemSetup, class tRunData > static __device__
			auto enact( typename tSystemSetup::Particle const&, tSystemSetup const&, tRunData const& aRunData ) -> typename tSystemSetup::value_type
			{
				return aRunData.randWalkStddev.value;
			}
		};
	}

	template< class tSystemSetup, class tRunData > __device__ inline
	auto rand_walk_stddev( typename tSystemSetup::Particle const& aParticle, tSystemSetup const& aSystem, tRunData const& aRunData ) -> typename tSystemSetup::value_type
	{
		return detail::RandWalkStddevImpl_<
			typename tSystemSetup::CompCount
		>::enact( aParticle, aSystem, aRunData );
	}
}
