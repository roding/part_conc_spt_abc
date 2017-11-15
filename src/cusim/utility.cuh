#ifndef UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6
#define UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6

#include "dsinline.cuh"

namespace cusim
{
	template< typename tReal, class tCompCount > __device__
	unsigned wrand_index( tReal aValue, tCompCount, DSInline<tReal,tCompCount> const& );


	template< class tSystemSetup, class tRunData, class tRand, class tRandData > __device__
	auto make_particle( tSystemSetup const&, tRunData const&, unsigned, tRand&, tRandData& )
		-> typename tSystemSetup::Particle;

	template< class tSystemSetup, class tRunData > __device__
	bool detected( typename tSystemSetup::Particle const&, tSystemSetup const&, tRunData const& );

	template< class tSystemSetup, class tRunData > __device__
	auto rand_walk_stddev( typename tSystemSetup::Particle const&, tSystemSetup const&, tRunData const& )
		-> typename tSystemSetup::value_type;
}

#include "utility.inl"
#endif // UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6
