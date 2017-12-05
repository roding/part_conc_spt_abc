#ifndef UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6
#define UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6

#include "dsinline.cuh"

namespace cusim
{
	template< typename tReal, class tCompCount > __device__
	unsigned wrand_index( tReal aValue, tCompCount, DSInline<tReal,tCompCount> const& );


	template< class tSimSetup, class tSimRun, class tRand, class tRandData > __device__
	auto make_particle( tSimSetup const&, tSimRun const&, unsigned, tRand&, tRandData& )
		-> typename tSimSetup::Particle;

	template< class tSimSetup, class tSimRun > __device__
	bool detected( typename tSimSetup::Particle const&, tSimSetup const&, tSimRun const& );

	template< class tSimSetup, class tSimRun, class tPart = typename tSimSetup::Particle > __device__
	auto rand_walk_stddev( tPart const&, tSimSetup const&, tSimRun const& )
		-> typename tSimSetup::Scalar;
}

#include "utility.inl"
#endif // UTILITY_CUH_5BBB274B_06CD_41D3_9BB9_5DA3E76827F6
