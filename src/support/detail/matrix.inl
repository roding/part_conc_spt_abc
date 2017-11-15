#include <algorithm>

#include <cassert>

#include "../strided_iterator.hpp"

namespace aspect
{
	struct MatrixRowMajor
	{
		template< typename tIter, class tNSize, class tMSize > 
		using ColIterator = StridedIterator<tIter, tMSize>;
		template< typename tIter, class, class >
		using RowIterator = tIter;
		
		template< class tMatrix > static constexpr
		std::size_t to_linear( tMatrix const& aMat, std::size_t aI, std::size_t aJ )
		{
			return aI * aMat.m() + aJ;
		}

		template< class tView, class tMatrix > static constexpr
		tView col_view( tMatrix&& aMat, std::size_t aCol )
		{
			using It_ = typename tView::iterator;
			return tView(
				It_( aMat.data() + aCol, aMat.m() ),
				It_( aMat.data() + aCol + aMat.n()*aMat.m(), aMat.m() )
			);
		}
		template< class tView, class tMatrix > static constexpr
		tView row_view( tMatrix&& aMat, std::size_t aRow )
		{
			using It_ = typename tView::iterator;
			return tView(
				It_( aMat.data() + aRow*aMat.m() ),
				It_( aMat.data() + aRow*aMat.m() + aMat.m() )
			);
		}
	};

	struct MatrixColMajor
	{
		template< typename tIter, class, class >
		using ColIterator = tIter;
		template< typename tIter, class tNSize, class tMSize > 
		using RowIterator = StridedIterator<tIter,tNSize>;

		template< class tMatrix > static constexpr
		std::size_t to_linear( tMatrix const& aMat, std::size_t aI, std::size_t aJ )
		{
			return aI + aMat.n() * aJ;
		}

		template< class tView, class tMatrix > static constexpr
		tView col_view( tMatrix&& aMat, std::size_t aCol )
		{
			using It_ = typename tView::iterator;
			return tView(
				It_( aMat.data() + aCol*aMat.n() ),
				It_( aMat.data() + aCol*aMat.n() + aMat.n() )
			);
		}
		template< class tView, class tMatrix > static constexpr
		tView row_view( tMatrix&& aMat, std::size_t aRow )
		{
			using It_ = typename tView::iterator;
			return tView(
				It_( aMat.data() + aRow, aMat.n() ),
				It_( aMat.data() + aRow + aMat.m()*aMat.n(), aMat.n() )
			);
		}
	};
}

namespace detail
{
	template< typename tIter, typename tConstIter > inline constexpr
	MatElementView<tIter,tConstIter>::MatElementView( tIter aBeg, tIter aEnd )
		: mBeg(aBeg)
		, mEnd(aEnd)
	{}

	template< typename tIter, typename tConstIter > inline CONSTEXPR_14
	auto MatElementView<tIter,tConstIter>::operator[] (std::size_t aIdx) -> value_type&
	{
		return mBeg[aIdx];
	}
	template< typename tIter, typename tConstIter > inline constexpr
	auto MatElementView<tIter,tConstIter>::operator[] (std::size_t aIdx) const -> value_type const&
	{
		return mBeg[aIdx];
	}

	template< typename tIter, typename tConstIter > inline constexpr
	bool MatElementView<tIter,tConstIter>::empty() const
	{
		return mEnd == mBeg;
	}
	template< typename tIter, typename tConstIter > inline constexpr
	std::size_t MatElementView<tIter,tConstIter>::size() const
	{
		return mEnd-mBeg;
	}

	template< typename tIter, typename tConstIter > inline CONSTEXPR_14
	tIter MatElementView<tIter,tConstIter>::begin()
	{
		return mBeg;
	}
	template< typename tIter, typename tConstIter > inline CONSTEXPR_14
	tIter MatElementView<tIter,tConstIter>::end()
	{
		return mEnd;
	}
	template< typename tIter, typename tConstIter > inline constexpr
	tConstIter MatElementView<tIter,tConstIter>::begin() const
	{
		return mBeg;
	}
	template< typename tIter, typename tConstIter > inline constexpr
	tConstIter MatElementView<tIter,tConstIter>::end() const
	{
		return mEnd;
	}
}

template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
Matrix<tScalar,tLayout,tN,tM>::Matrix()
	: mData( nullptr )
	, mN()
	, mM()
{}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
Matrix<tScalar,tLayout,tN,tM>::Matrix( tScalar* aData, NSizeType aN, MSizeType aM )
	: mData( aData )
	, mN( aN )
	, mM( aM )
{}

template< typename tScalar, class tLayout, class tN, class tM > inline
Matrix<tScalar,tLayout,tN,tM>::Matrix( NSizeType aN, MSizeType aM )
	: mData( new tScalar[aN*aM] )
	, mN( aN )
	, mM( aM )
{}

template< typename tScalar, class tLayout, class tN, class tM > inline
Matrix<tScalar,tLayout,tN,tM>::Matrix( NSizeType aN, MSizeType aM, tScalar aFill )
	: Matrix<tScalar,tLayout,tN,tM>( aN, aM )
{
	std::fill_n( mData, mN*mM, aFill );
}


template< typename tScalar, class tLayout, class tN, class tM > inline
Matrix<tScalar,tLayout,tN,tM>::~Matrix()
{
	delete [] mData;
}


template< typename tScalar, class tLayout, class tN, class tM > inline
Matrix<tScalar,tLayout,tN,tM>::Matrix( Matrix const& aOther )
	: Matrix<tScalar,tLayout,tN,tM>( aOther.n(), aOther.m() )
{
	std::copy_n( mData, mN*mM, aOther.data() );
}
template< typename tScalar, class tLayout, class tN, class tM > inline
Matrix<tScalar,tLayout,tN,tM>& Matrix<tScalar,tLayout,tN,tM>::operator= (Matrix const& aOther)
{
	delete [] mData;

	mN = aOther.n();
	mM = aOther.m();

	mData = new tScalar[aOther.n()*aOther.m()];
	std::copy_n( mData, mN*mM, aOther.data() );
}


template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
tScalar& Matrix<tScalar,tLayout,tN,tM>::operator[] (std::size_t aIdx)
{
	return assert( aIdx < mN*mM ), mData[aIdx];
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
tScalar const& Matrix<tScalar,tLayout,tN,tM>::operator[] (std::size_t aIdx) const
{
	return assert( aIdx < mN*mM ), mData[aIdx];
}

template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
tScalar& Matrix<tScalar,tLayout,tN,tM>::operator() (std::size_t aI, std::size_t aJ)
{
	return assert( aI < mN && aJ < mM ), mData[ tLayout::to_linear(*this,aI,aJ) ];
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
tScalar const& Matrix<tScalar,tLayout,tN,tM>::operator() (std::size_t aI, std::size_t aJ) const
{
	return assert( aI < mN && aJ < mM ), mData[ tLayout::to_linear(*this,aI,aJ) ];
}


template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
tScalar* Matrix<tScalar,tLayout,tN,tM>::data()
{
	return mData;
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
tScalar const* Matrix<tScalar,tLayout,tN,tM>::data() const
{
	return mData;
}

template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::n() const -> NSizeType
{
	return mN;
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::m() const -> MSizeType
{
	return mM;
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::size() const -> SizeType
{
	return mN*mM;
}

template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
auto Matrix<tScalar,tLayout,tN,tM>::row( std::size_t aRow ) -> RowView
{
	return Layout::template row_view<RowView>( *this, aRow );
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::row( std::size_t aRow ) const -> ConstRowView
{
	return Layout::template row_view<ConstRowView>( *this, aRow );
}

template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
auto Matrix<tScalar,tLayout,tN,tM>::col( std::size_t aCol ) -> ColView
{
	return Layout::template col_view<ColView>( *this, aCol );
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::col( std::size_t aCol ) const -> ConstColView
{
	return Layout::template col_view<ConstColView>( *this, aCol );
}

template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
auto Matrix<tScalar,tLayout,tN,tM>::linear_elements() -> ElementView
{
	return ElementView( mData, mData+mN*mM );
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::linear_elements() const -> ConstElementView
{
	return ConstElementView( mData, mData+mN*mM );
}


template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
auto Matrix<tScalar,tLayout,tN,tM>::lend() -> tScalar*
{
	return mData+mN*mM;
}
template< typename tScalar, class tLayout, class tN, class tM > inline CONSTEXPR_14
auto Matrix<tScalar,tLayout,tN,tM>::lbegin() -> tScalar*
{
	return mData;
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::lend() const -> tScalar const*
{
	return mData+mN*mM;
}
template< typename tScalar, class tLayout, class tN, class tM > inline constexpr
auto Matrix<tScalar,tLayout,tN,tM>::lbegin() const -> tScalar const*
{
	return mData;
}
