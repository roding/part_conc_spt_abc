#include <random>
#include <algorithm>
#include <functional> // std::ref
#include <type_traits>

template< class tRng, std::size_t tStateSize > inline
tRng make_prng()
{
	std::random_device rsrc{};
	
	typename tRng::result_type data[tStateSize];
	std::generate( data, data+tStateSize, std::ref(rsrc) );

	std::seed_seq sseq( data, data+tStateSize );
	return tRng( sseq );
}

template< typename tTo, typename tFrom > inline
tTo int_cast( tFrom aX )
{
	auto const y = static_cast<tTo>(aX);

	if( static_cast<tFrom>(y) != aX )
		detail::throw_bad_int_cast();

	constexpr bool same_ = std::is_signed<tTo>::value == std::is_signed<tFrom>::value;
	if( !same_ && (y < tTo{}) != (aX < tFrom{}) )
		detail::throw_bad_int_cast();

	return y;
}
