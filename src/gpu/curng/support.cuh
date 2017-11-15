/*-******************************************************* -- HEADER -{{{1- */
/*-	Support
 */
/*-***************************************************************** -}}}1- */

#ifndef SUPPORT_CUH_D3250AEF_9C00_48C3_B55F_C873629BF67C
#define SUPPORT_CUH_D3250AEF_9C00_48C3_B55F_C873629BF67C

#include <cstddef>

namespace curng
{
	/** Generate uniform real ∈ [0,1)
	 *
	 * Generate uniform real value \f$x \in \left[0, 1\right)\f$. This code is
     * roughly based on the implementation in libc++, available at
	 *
	 *   https://llvm.org/svn/llvm-project/libcxx/trunk/include/random
	 *
	 * \note The above code isn't strictly conforming to the standard at the
	 * time of writing (revision r31598) as the implementation produces numbers
	 * in the closed range \f$\left[0,1\right]\f$. GCC's libstdc++ has a
	 * similar problem, and fixes this by rejecting any numbers \f$\ge 1\f$ and
	 * retrying. This includes a fix in the same spirit. There is a small but
	 * measurable overhead (the rejection method leads to occasional warp
	 * divergence).
	 */
	template< typename tReal, std::size_t tBits, class tEngine, class tEngData >
	__device__
	tReal generate_canonical( tEngine&, unsigned, tEngData& );
}

#include "support.inl"
#endif // SUPPORT_CUH_D3250AEF_9C00_48C3_B55F_C873629BF67C
