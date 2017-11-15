/*-******************************************************* -- HEADER -{{{1- */
/*-	Strided iterator
 */
/*-***************************************************************** -}}}1- */

#ifndef STRIDED_ITERATOR_HPP_78A4CA4C_FA02_41E8_B807_388C67B9F4BB
#define STRIDED_ITERATOR_HPP_78A4CA4C_FA02_41E8_B807_388C67B9F4BB

#include <iterator>
#include <type_traits>

#include <cstddef>

#include "compat.hpp"
#include "dsvalue.hpp"

template< typename tIter, class tStrideValue >
class StridedIterator
{
	static_assert( std::is_same<typename std::iterator_traits<tIter>::iterator_category, std::random_access_iterator_tag>::value, "StridedIterator: requires RandomAccessIterator" );
	
	public:
		explicit constexpr 
		StridedIterator( tIter, tStrideValue = tStrideValue() );

		StridedIterator() = default;

		StridedIterator( StridedIterator const& ) = default;
		StridedIterator& operator= (StridedIterator const&) = default;

		StridedIterator( StridedIterator&& ) = default;
		StridedIterator& operator= (StridedIterator&&) = default;

	public:
		using difference_type = typename std::iterator_traits<tIter>::difference_type;
		using value_type = typename std::iterator_traits<tIter>::value_type;
		using pointer = typename std::iterator_traits<tIter>::pointer;
		using reference = typename std::iterator_traits<tIter>::reference;
		using iterator_category = std::random_access_iterator_tag;

	public:
		constexpr pointer operator-> () const;
		constexpr reference operator* () const;
		constexpr reference operator[] (difference_type) const;

		CONSTEXPR_14 StridedIterator& operator++();
		constexpr StridedIterator operator++(int) const;
		CONSTEXPR_14 StridedIterator& operator--();
		constexpr StridedIterator operator--(int) const;

		constexpr StridedIterator operator+ (difference_type) const;
		constexpr StridedIterator operator- (difference_type) const;
		CONSTEXPR_14 StridedIterator& operator+= (difference_type);
		CONSTEXPR_14 StridedIterator& operator-= (difference_type);

		constexpr difference_type operator- (StridedIterator const&) const;

		constexpr bool operator== (StridedIterator const& aOther) const;
		constexpr bool operator!= (StridedIterator const& aOther) const;
		constexpr bool operator< (StridedIterator const& aOther) const;
		constexpr bool operator> (StridedIterator const& aOther) const;
		constexpr bool operator<= (StridedIterator const& aOther) const;
		constexpr bool operator>= (StridedIterator const& aOther) const;
		
	private:
		tIter mIt;
		tStrideValue mStride;
};

template< typename tIter, class tStride > constexpr
auto operator+ (std::ptrdiff_t, StridedIterator<tIter,tStride> const&)
	-> StridedIterator<tIter,tStride>;
template< typename tIter, class tStride > constexpr
auto operator- (std::ptrdiff_t, StridedIterator<tIter,tStride> const&)
	-> StridedIterator<tIter,tStride>;

#include "detail/strided_iterator.inl"
#endif // STRIDED_ITERATOR_HPP_78A4CA4C_FA02_41E8_B807_388C67B9F4BB
