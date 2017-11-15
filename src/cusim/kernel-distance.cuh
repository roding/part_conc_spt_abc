#ifndef KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161
#define KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161

namespace cusim
{
	// N and M are expected to be relatively small
	
	template< 
		typename tReal,
		typename tCount
	>
	__global__
	void K_distance( 
		tReal* aDistance,
		unsigned aN, 
		unsigned aM, 
		tCount* aCur,  // WARNING: destructive (in-place)
		tCount const* aSATReference
	);
}

#include "kernel-distance.inl"
#endif // KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161
