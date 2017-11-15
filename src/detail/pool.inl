#include "support/cuda_error.hpp"

#include <cassert>

template< typename tType > inline
Pool<tType>::Pool( std::size_t aArrayCount, std::size_t aInitial )
	: mFreePool( aInitial, nullptr )
	, mArrayCount(aArrayCount)
	, mAllocated(0)
{
	for( std::size_t i = 0; i < mFreePool.size(); ++i )
		CUDA_CHECKED cudaMalloc( &mFreePool[i], sizeof(tType)*mArrayCount );
}

template< typename tType > inline
Pool<tType>::~Pool() noexcept
{
	if( mAllocated )
	{
		std::fprintf( stderr, "Warning: Pool: %zu allocations remain.\n", mAllocated );
	}

	for( auto ptr : mFreePool )
		cudaFree( ptr );
}


template< typename tType > inline
Pool<tType>::Pool( Pool&& aOther ) noexcept
	: mFreePool( std::move(aOther.mFreePool) )
	, mArrayCount(aOther.mArrayCount)
	, mAllocated(aOther.mAllocated)
{
	aOther.mFreePool.clear();
	aOther.mAllocated = 0;
}
template< typename tType > inline
auto Pool<tType>::operator= (Pool&& aOther) noexcept -> Pool&
{
	std::swap( mFreePool, aOther.mFreePool );
	std::swap( mArrayCount, aOther.mArrayCount );
	std::swap( mAllocated, aOther.mAllocated );
	return *this;
}


template< typename tType > inline
tType* Pool<tType>::alloc()
{
	if( !mFreePool.empty() )
	{
		auto ret = mFreePool.back();
		mFreePool.pop_back();
		++mAllocated;
		return ret;
	}

	// linear allocation policy. :-(
	tType* ptr = nullptr;
	CUDA_CHECKED cudaMalloc( &ptr, sizeof(tType)*mArrayCount );

	++mAllocated;
	return ptr;
}
template< typename tType > inline
void Pool<tType>::free( tType const* aPtr )
{
	assert( aPtr );
	assert( mAllocated > 0 );

	mFreePool.push_back( const_cast<tType*>(aPtr) );
	--mAllocated;
}

template< typename tType > inline
void Pool<tType>::free_cuda_resources()
{
	for( auto ptr : mFreePool )
		cudaFree( ptr );

	mFreePool.clear();
}
