#include "../support/utility.hpp"
#include "../support/cuda_error.hpp"

namespace cusim
{
	template< typename tCount > __host__ inline
	Histogram2D<tCount>::Histogram2D( std::size_t aN, std::size_t aM, Count const* aPtr )
		: mN(aN)
		, mM(aM)
		, mDevBuffer(nullptr)
	{
		CUDA_CHECKED cudaMalloc( &mDevBuffer, aN*aM*sizeof(Count) );

		if( aPtr )
		{
			CUDA_CHECKED cudaMemcpy( mDevBuffer, aPtr, aN*aM*sizeof(Count), cudaMemcpyHostToDevice );
		}
	}

	template< typename tCount > __host__ inline
	Histogram2D<tCount>::~Histogram2D()
	{
		if( mDevBuffer );
			cudaFree( mDevBuffer );
	}


	template< typename tCount > __host__ inline
	Histogram2D<tCount>::Histogram2D( Histogram2D&& aOther ) noexcept
		: mN(aOther.mN)
		, mM(aOther.mM)
		, mDevBuffer(aOther.mDevBuffer)
	{
		aOther.mN = aOther.mM = 0;
		aOther.mDevBuffer = nullptr;
	}
	template< typename tCount > __host__ inline
	auto Histogram2D<tCount>::operator= (Histogram2D&& aOther) noexcept -> Histogram2D&
	{
		std::swap( mN, aOther.mN );
		std::swap( mM, aOther.mM );
		std::swap( mDevBuffer, aOther.mDevBuffer );
		return *this;
	}

	template< typename tCount > __host__ inline
	void Histogram2D<tCount>::clear()
	{
		CUDA_CHECKED cudaMemset( mDevBuffer, 0, mN*mM*sizeof(Count) );
	}
	template< typename tCount > __host__ inline
	void Histogram2D<tCount>::clear_async( cudaStream_t aStream )
	{
		CUDA_CHECKED cudaMemsetAsync( mDevBuffer, 0, mN*mM*sizeof(Count), aStream );
	}

	template< typename tCount > __host__ inline
	auto Histogram2D<tCount>::device_ptr() -> Count*
	{
		return mDevBuffer;
	}
	template< typename tCount > __host__ inline
	auto Histogram2D<tCount>::device_ptr() const -> Count const*
	{
		return mDevBuffer;
	}

	template< typename tCount > __host__ inline
	std::size_t Histogram2D<tCount>::n() const
	{
		return mN;
	}
	template< typename tCount > __host__ inline
	std::size_t Histogram2D<tCount>::m() const
	{
		return mM;
	}

	template< typename tCount > __host__ inline
	void Histogram2D<tCount>::free_cuda_resources()
	{
		if( mDevBuffer ) 
		{
			cudaFree( mDevBuffer );
			mDevBuffer = nullptr;
		}
	}



	template< typename tCount, typename tReal > __host__ inline
	HistogramRecorder<tCount,tReal>::HistogramRecorder( Histogram2D<tCount>& aHist, tReal aDDE )
		: mDDE(aDDE)
		, mDEBinCount(int_cast<tCount>(aHist.m()))
		, mHistoBuffer(aHist.device_ptr())
	{}

	template< typename tCount, typename tReal > __device__ inline
	void HistogramRecorder<tCount,tReal>::record( tCount aK, tReal aDE )
	{
		auto const de = tCount(aDE / mDDE);
		if( de < mDEBinCount )
		{
			atomicInc( &mHistoBuffer[aK*mDEBinCount+de], std::numeric_limits<tCount>::max() );
		}
	}
}
