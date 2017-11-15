#include "support/cuda_error.hpp"

#include <cassert>

inline
EvPool::EvPool( std::size_t aInitial )
	: mFreePool( aInitial, nullptr )
	, mAllocated(0)
{
	for( std::size_t i = 0; i < mFreePool.size(); ++i )
		CUDA_CHECKED cudaEventCreate( &mFreePool[i] );
}

inline
EvPool::~EvPool() noexcept
{
	if( mAllocated )
	{
		std::fprintf( stderr, "Warning: EvPool: %zu allocations remain.\n", mAllocated );
	}

	for( auto ev : mFreePool )
		cudaEventDestroy( ev );
}


inline
EvPool::EvPool( EvPool&& aOther ) noexcept
	: mFreePool( std::move(aOther.mFreePool) )
	, mAllocated(aOther.mAllocated)
{
	aOther.mFreePool.clear();
	aOther.mAllocated = 0;
}
inline
auto EvPool::operator= (EvPool&& aOther) noexcept -> EvPool&
{
	std::swap( mFreePool, aOther.mFreePool );
	std::swap( mAllocated, aOther.mAllocated );
	return *this;
}


inline
cudaEvent_t EvPool::alloc()
{
	if( !mFreePool.empty() )
	{
		auto ret = mFreePool.back();
		mFreePool.pop_back();
		++mAllocated;
		return ret;
	}

	// linear allocation policy. :-(
	cudaEvent_t ev;
	CUDA_CHECKED cudaEventCreate( &ev );

	++mAllocated;
	return ev;
}
inline
void EvPool::free( cudaEvent_t aEvent )
{
	assert( mAllocated > 0 );
	mFreePool.push_back( aEvent );
	--mAllocated;
}

inline
void EvPool::free_cuda_resources()
{
	for( auto ev : mFreePool )
		cudaEventDestroy( ev );

	mFreePool.clear();
}
