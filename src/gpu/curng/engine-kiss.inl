#include <limits>
#include <vector>
#include <random>

#include "../support/cuda_error.hpp"

namespace curng
{
	__device__ inline
	EngineKISS::Engine::Engine( unsigned aTid, GlobalData const& aData )
		: jsr(aData.jsr[aTid])
		, z(aData.z[aTid])
		, w(aData.w[aTid])
		, jcong(aData.jcong[aTid])
	{}

	__device__ inline
	auto EngineKISS::Engine::operator() (unsigned, GlobalData&) -> result_type
	{
		auto const jz = jsr;
		jsr ^= jsr << 13;
		jsr ^= jsr >> 17;
		jsr ^= jsr << 5;

		z = 36969 * (z & 65535) + (z>>16);
		w = 18000 * (w & 65535) + (w>>16);

		jcong = 69069 * jcong + 1234567;

		return (((z << 16) + w) ^ jcong) + (jz+jsr);
	}

	__device__ inline
	void EngineKISS::Engine::store( unsigned aTid, GlobalData& aData )
	{
		aData.jsr[aTid] = jsr;
		aData.z[aTid] = z;
		aData.w[aTid] = w;
		aData.jcong[aTid] = jcong;
	}

	__host__ __device__ constexpr
	auto EngineKISS::Engine::min() -> result_type
	{
		return std::numeric_limits<result_type>::min();
	}
	__host__ __device__ constexpr
	auto EngineKISS::Engine::max() -> result_type
	{
		return std::numeric_limits<result_type>::max();
	}

	template< class tHostRNG > __host__ inline
	auto EngineKISS::initialize( std::size_t aThreadCount, tHostRNG& aRng ) -> GlobalData
	{
		std::vector<std::uint32_t> buffer( aThreadCount );
		std::uniform_int_distribution<std::uint32_t> dist( 0, std::numeric_limits<std::uint32_t>::max() );

		std::size_t const bytes = sizeof(std::uint32_t)*aThreadCount;
		
		EngineKISS::GlobalData ret;
		CUDA_CHECKED cudaMalloc( &ret.jsr, bytes );
		std::generate( buffer.begin(), buffer.end(), [&] { return dist(aRng); } );
		CUDA_CHECKED cudaMemcpy( ret.jsr, buffer.data(), bytes, cudaMemcpyHostToDevice );

		CUDA_CHECKED cudaMalloc( &ret.z, bytes );
		std::generate( buffer.begin(), buffer.end(), [&] { return dist(aRng); } );
		CUDA_CHECKED cudaMemcpy( ret.z, buffer.data(), bytes, cudaMemcpyHostToDevice );

		CUDA_CHECKED cudaMalloc( &ret.w, bytes );
		std::generate( buffer.begin(), buffer.end(), [&] { return dist(aRng); } );
		CUDA_CHECKED cudaMemcpy( ret.w, buffer.data(), bytes, cudaMemcpyHostToDevice );

		CUDA_CHECKED cudaMalloc( &ret.jcong, bytes );
		std::generate( buffer.begin(), buffer.end(), [&] { return dist(aRng); } );
		CUDA_CHECKED cudaMemcpy( ret.jcong, buffer.data(), bytes, cudaMemcpyHostToDevice );

		return ret;
	}

	__host__ inline
	void EngineKISS::cleanup( GlobalData& aData )
	{
		if( aData.jsr ) cudaFree( aData.jsr );
		if( aData.z ) cudaFree( aData.z );
		if( aData.w ) cudaFree( aData.w );
		if( aData.jcong ) cudaFree( aData.jcong );

		aData.jsr = aData.z = aData.w = aData.jcong = nullptr;
	}
}
