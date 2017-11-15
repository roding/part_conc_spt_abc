#include "../support/utility.hpp"
#include "../support/cuda_error.hpp"

namespace cusim
{
	template< typename tCount, typename tReal > __device__ inline
	void OutputHisto2D<tCount,tReal>::Recorder::record( tCount aK, tReal aDE )
	{
		tCount const de = unsigned(aDE / dDE);
		if( de < deBinCount )
			atomicInc( &buffer[aK*deBinCount + de], std::numeric_limits<tCount>::max() );
	}
	

	template< typename tCount, typename tReal > __host__ inline
	OutputHisto2D<tCount,tReal>::OutputHisto2D( std::size_t aKs, std::size_t aDEs, tReal aDDE, tCount const* aInitial )
		: mKs(aKs)
		, mDEs(aDEs)
		, mDDE(aDDE)
		, mDevBuffer( nullptr )
	{
		CUDA_CHECKED cudaMalloc( &mDevBuffer, sizeof(tCount)*mKs*mDEs );

		if( aInitial )
		{
			CUDA_CHECKED cudaMemcpy( mDevBuffer, aInitial, sizeof(tCount)*mKs*mDEs, cudaMemcpyHostToDevice );
		}
	}

	template< typename tCount, typename tReal > __host__ inline
	OutputHisto2D<tCount,tReal>::~OutputHisto2D()
	{
		cudaFree( mDevBuffer );
	}


	template< typename tCount, typename tReal > __host__ inline
	void OutputHisto2D<tCount,tReal>::reset()
	{
		CUDA_CHECKED cudaMemset( mDevBuffer, 0, sizeof(tCount)*mKs*mDEs );
	}
	template< typename tCount, typename tReal > __host__ inline
	auto OutputHisto2D<tCount,tReal>::recorder() const -> Recorder
	{
		return Recorder{
			int_cast<tCount>(mDEs),
			mDDE,
			mDevBuffer
		};
	}
}
