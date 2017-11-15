#include "../support/cuda_error.hpp"

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
	if( mAllocated )
	{
		std::fprintf( stderr, "Warning: Pool::free_cuda_resources(): %zu allocations remain.\n", mAllocated );
	}

	for( auto ptr : mFreePool )
		cudaFree( ptr );

	mFreePool.clear();
}



template< typename tType > inline
MappedPool<tType>::MappedPool( std::size_t aArrayCount, std::size_t aInitial )
	: mFreePoolHost( aInitial, nullptr )
	, mFreePoolDev( aInitial, nullptr )
	, mArrayCount(aArrayCount)
	, mAllocated(0)
{
	for( std::size_t i = 0; i < mFreePoolHost.size(); ++i )
	{
		CUDA_CHECKED cudaHostAlloc( &mFreePoolHost[i], sizeof(tType)*mArrayCount, cudaHostAllocMapped );
		CUDA_CHECKED cudaHostGetDevicePointer( &mFreePoolDev[i], mFreePoolHost[i], 0 );
	}
}

template< typename tType > inline
MappedPool<tType>::~MappedPool() noexcept
{
	if( mAllocated )
	{
		std::fprintf( stderr, "Warning: MappedPool: %zu allocations remain.\n", mAllocated );
	}

	for( auto ptr : mFreePoolHost )
		cudaFreeHost( ptr );
}


template< typename tType > inline
MappedPool<tType>::MappedPool( MappedPool&& aOther ) noexcept
	: mFreePoolHost( std::move(aOther.mFreePoolHost) )
	, mFreePoolDev( std::move(aOther.mFreePoolDev) )
	, mArrayCount(aOther.mArrayCount)
	, mAllocated(aOther.mAllocated)
{
	aOther.mFreePoolDev.clear();
	aOther.mFreePoolHost.clear();
	aOther.mAllocated = 0;
}
template< typename tType > inline
auto MappedPool<tType>::operator= (MappedPool&& aOther) noexcept -> MappedPool&
{
	std::swap( mFreePoolDev, aOther.mFreePoolDev );
	std::swap( mFreePoolHost, aOther.mFreePoolHost );
	std::swap( mArrayCount, aOther.mArrayCount );
	std::swap( mAllocated, aOther.mAllocated );
	return *this;
}


template< typename tType > inline
std::tuple<tType*,tType*> MappedPool<tType>::alloc()
{
	if( !mFreePoolHost.empty() )
	{
		auto hret = mFreePoolHost.back();
		auto dret = mFreePoolDev.back();
		mFreePoolHost.pop_back();
		mFreePoolDev.pop_back();
		++mAllocated;
		return std::make_tuple(hret,dret);
	}

	// linear allocation policy. :-(
	tType* hptr = nullptr;
	CUDA_CHECKED cudaHostAlloc( &hptr, sizeof(tType)*mArrayCount, cudaHostAllocMapped );

	tType* dptr = nullptr;
	CUDA_CHECKED cudaHostGetDevicePointer( &dptr, hptr, 0 );

	++mAllocated;
	return std::make_tuple(hptr,dptr);
}
template< typename tType > inline
void MappedPool<tType>::free( tType const* aHost, tType const* aDev )
{
	assert( aHost && aDev );
	assert( mAllocated > 0 );

	mFreePoolHost.push_back( const_cast<tType*>(aHost) );
	mFreePoolDev.push_back( const_cast<tType*>(aDev) );
	--mAllocated;
}

template< typename tType > inline
void MappedPool<tType>::free_cuda_resources()
{
	if( mAllocated )
	{
		std::fprintf( stderr, "Warning: MappedPool::free_cuda_resources(): %zu allocations remain.\n", mAllocated );
	}

	for( auto ptr : mFreePoolHost )
		cudaFreeHost( ptr );

	mFreePoolHost.clear();
	mFreePoolDev.clear();
}
