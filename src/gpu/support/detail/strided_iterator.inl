template< typename tIter, class tStride > inline constexpr
StridedIterator<tIter,tStride>::StridedIterator( tIter aBase, tStride aStride )
	: mIt(aBase)
	, mStride(aStride)
{}

template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator->() const -> pointer
{
	return &*mIt;
}
template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator*() const -> reference
{
	return *mIt;
}
template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator[] (difference_type aIdx) const -> reference
{
	return mIt[aIdx * difference_type(mStride)];
}

template< typename tIter, class tStride > inline CONSTEXPR_14
auto StridedIterator<tIter,tStride>::operator++() -> StridedIterator&
{
	mIt += mStride;
	return *this;
}
template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator++(int) const -> StridedIterator
{
	return StridedIterator{ mIt + mStride, mStride };
}
template< typename tIter, class tStride > inline CONSTEXPR_14
auto StridedIterator<tIter,tStride>::operator--() -> StridedIterator&
{
	mIt -= mStride;
	return *this;
}
template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator--(int) const -> StridedIterator
{
	return StridedIterator{ mIt - mStride, mStride };
}

template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator+ (difference_type aX) const -> StridedIterator
{
	return StridedIterator( mIt + aX*difference_type(mStride), mStride );
}
template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator- (difference_type aX) const -> StridedIterator
{
	return StridedIterator( mIt - aX*difference_type(mStride), mStride );
}

template< typename tIter, class tStride > inline CONSTEXPR_14
auto StridedIterator<tIter,tStride>::operator+= (difference_type aX) -> StridedIterator&
{
	mIt += aX * difference_type(mStride);
	return *this;
}
template< typename tIter, class tStride > inline CONSTEXPR_14
auto StridedIterator<tIter,tStride>::operator-= (difference_type aX) -> StridedIterator&
{
	mIt -= aX * difference_type(mStride);
	return *this;
}

template< typename tIter, class tStride > inline constexpr
auto StridedIterator<tIter,tStride>::operator- (StridedIterator const& aOther) const -> difference_type
{
	return mIt - aOther.mIt;
}

template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator== (StridedIterator const& aOther) const
{
	return mIt == aOther.mIt; // NOTE: we're ignoring the stride here!
}
template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator!= (StridedIterator const& aOther) const
{
	return mIt != aOther.mIt;
}
template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator< (StridedIterator const& aOther) const
{
	return mIt < aOther.mIt;
}
template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator> (StridedIterator const& aOther) const
{
	return mIt > aOther.mIt;
}
template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator<= (StridedIterator const& aOther) const
{
	return mIt <= aOther.mIt;
}
template< typename tIter, class tStride > inline constexpr
bool StridedIterator<tIter,tStride>::operator>= (StridedIterator const& aOther) const
{
	return mIt >= aOther.mIt;
}

template< typename tIter, class tStride > constexpr
auto operator+ (std::ptrdiff_t aX, StridedIterator<tIter,tStride> const& aIt)
	-> StridedIterator<tIter,tStride>
{
	return aIt + aX;
}
template< typename tIter, class tStride > constexpr
auto operator- (std::ptrdiff_t aX, StridedIterator<tIter,tStride> const& aIt)
	-> StridedIterator<tIter,tStride> 
{
	return aIt - aX;
}

