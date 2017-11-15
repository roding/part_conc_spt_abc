#include <limits>
#include <vector>
#include <random>

#include "../support/cuda_error.hpp"

namespace curng
{
	__device__ inline
	EngineLCG48_32::Engine::Engine( unsigned aTid, GlobalData const& aData )
		: mState(aData.states[aTid])
	{}

	__device__ inline
	auto EngineLCG48_32::Engine::operator() (unsigned, GlobalData&) -> result_type
	{
		static constexpr auto kMask_ = (1ull<<48)-1ull;
		mState = (a * mState + c) & kMask_;
		return result_type(mState >> 16);
	}

	__device__ inline
	void EngineLCG48_32::Engine::store( unsigned aTid, GlobalData& aData )
	{
		aData.states[aTid] = mState;
	}

	__host__ __device__ constexpr
	auto EngineLCG48_32::Engine::min() -> result_type
	{
		return std::numeric_limits<result_type>::min();
	}
	__host__ __device__ constexpr
	auto EngineLCG48_32::Engine::max() -> result_type
	{
		return std::numeric_limits<result_type>::max();
	}

	template< class tHostRNG > __host__ inline
	auto EngineLCG48_32::initialize( std::size_t aThreadCount, tHostRNG& aRng ) -> GlobalData
	{
		std::vector<std::uint64_t> buffer( aThreadCount );
		std::uniform_int_distribution<std::uint64_t> dist( 0, (1ull<<48)-1 );

		std::size_t const bytes = sizeof(std::uint64_t)*aThreadCount;
		
		EngineLCG48_32::GlobalData ret;
		CUDA_CHECKED cudaMalloc( &ret.states, bytes );
		std::generate( buffer.begin(), buffer.end(), [&] { return dist(aRng); } );
		CUDA_CHECKED cudaMemcpy( ret.states, buffer.data(), bytes, cudaMemcpyHostToDevice );
		return ret;
	}

	__host__ inline
	void EngineLCG48_32::cleanup( GlobalData& aData )
	{
		if( aData.states ) cudaFree( aData.states );

		aData.states = nullptr;
	}
}
