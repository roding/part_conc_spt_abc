/*-******************************************************* -- HEADER -{{{1- */
/*-	Static / Dynamic "vector"
 *
 * Pick either std::vector (dynamic) or std::array (static).
 */
/*-***************************************************************** -}}}1- */

#ifndef DSVECTOR_HPP_49B95419_A46D_477F_BC63_367FFAE9C587
#define DSVECTOR_HPP_49B95419_A46D_477F_BC63_367FFAE9C587

#include <array>
#include <vector>

#include "dsvalue.hpp"

namespace detail
{
	template< class > 
	struct DVector
	{
		template< typename tType >
		using container = std::vector<tType>;
	};
	template< class tSVal > 
	struct SVector
	{
		static constexpr typename tSVal::type size_ = tSVal::staticValue;
		template< typename tType >
		using container = std::array<tType,size_>;
	};
}

template< typename tElementType, class tSizeValue >
using DSVector = typename std::conditional<
	is_static_value<tSizeValue>::value,
	detail::SVector<tSizeValue>,
	detail::DVector<tSizeValue>
>::type::template container<tElementType>;


//TODO: make_dsvector that does the initial resize (or uses a constructor?)


template< typename tEType, class tDSValue > inline
void resize( std::vector<tEType>& aVector, tDSValue const& aSize )
{
	aVector.resize( aSize );
}
template< typename tEType, std::size_t tASize, class tDSValue > inline
void resize( std::array<tEType,tASize>& aArray, tDSValue const& aSize )
{
	assert( aArray.size() == aSize );
}

template< typename tEType, class tDSValue > inline
void resize( std::vector<tEType>& aVector, tDSValue const& aSize, tEType const& aX )
{
	aVector.resize( aSize, aX );
}
template< typename tEType, std::size_t tASize, class tDSValue > inline
void resize( std::array<tEType,tASize>& aArray, tDSValue const& aSize, tEType const& )
{
	assert( aArray.size() == aSize );
}

#endif // DSVECTOR_HPP_49B95419_A46D_477F_BC63_367FFAE9C587
