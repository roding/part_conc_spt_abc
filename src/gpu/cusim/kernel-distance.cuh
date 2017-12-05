#ifndef KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161
#define KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161

namespace cusim
{
	/** Compute cumulative histogram distance
	 * 
	 * Compute the cumulative histogram distance.
	 *
	 * <b>IMPORTANT:</b>
	 *  - the reference \a aSATReference is a cumulative histogram
	 *  - the input histogram \a aCur will be reset to zero by `K_distance`
	 *
	 * `K_distance` assumes that the values of \a aN and \a aM are relatively
	 * small: a single block of \f$32\times32\f$ performs the computation. 32
	 * columns/rows are handled in parallel while 32 threads iterate over the
	 * elements of the column/row.
	 *
	 * \warning Assumes a warp size of 32 and assumes that the kernel is
	 * launched in a single 32xN block.
	 */
	template< 
		typename tReal,
		typename tCount
	>
	__global__
	void K_distance( 
		tReal* aDistance,
		unsigned aN, 
		unsigned aM, 
		tCount* aCur, // !NOTE! reset to zero by the kernel
		tCount const* aSATReference
	);

	extern template __global__
	void K_distance<float,unsigned>( float*, unsigned, unsigned, unsigned*, unsigned const* );
	extern template __global__
	void K_distance<double,unsigned>( double*, unsigned, unsigned, unsigned*, unsigned const* );
}

#include "kernel-distance.inl"
#endif // KERNEL_DISTANCE_CUH_4ED6C75E_9C8E_4226_882D_26B613CD3161
