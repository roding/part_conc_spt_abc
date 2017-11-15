/*-******************************************************* -- HEADER -{{{1- */
/*-	utility.hpp 
 */
/*-***************************************************************** -}}}1- */

#ifndef UTILITY_HPP_1609C312_1550_4C57_B557_B56AAA6A2DFC
#define UTILITY_HPP_1609C312_1550_4C57_B557_B56AAA6A2DFC

#include <cstddef>

template< typename tType >
struct Identity
{
	using type = tType;
};

/** Create seeded random engine
 *
 * Creates, seeds and returns an instance of the \a tRng random engine. Tries
 * to seed it properly (i.e., with more than a single `result_type` worth of
 * bits should the state be larger than this [1]), but probably fails to do so
 * completely unbiased [2].
 *
 * Uses `std::random_device` as a source for the seeds.
 *
 * [1] https://codereview.stackexchange.com/a/109266
 * [2] http://www.pcg-random.org/posts/cpp-seeding-surprises.html
 */
template< class tRng, std::size_t tStateSize = tRng::state_size >
tRng make_prng();

/** Checked cast between integer types
 *
 * C.f. `narrow()` of GSL. Converts from one integer type to another, and
 * throws an error if the cast changed the value.
 */
template< typename tTo, typename tFrom >
tTo int_cast( tFrom );

namespace detail
{
	[[noreturn]] void throw_bad_int_cast();
}

#include "detail/utility.inl"
#endif // UTILITY_HPP_1609C312_1550_4C57_B557_B56AAA6A2DFC
