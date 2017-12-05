#include "kernel-distance.cuh"

namespace cusim
{
	template __global__
	void K_distance<float,unsigned>( float*, unsigned, unsigned, unsigned*, unsigned const* );
	template __global__
	void K_distance<double,unsigned>( double*, unsigned, unsigned, unsigned*, unsigned const* );
}
