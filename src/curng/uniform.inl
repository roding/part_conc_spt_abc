#include <limits>

#include "support.cuh"

namespace curng
{
	template< typename tReal > __device__ inline
	UniformReal::Distribution<tReal>::Distribution( unsigned, UniformReal::GlobalData<tReal> const& )
	{}

	template< typename tReal > template< class tRand, class tRandData > __device__ inline
	auto UniformReal::Distribution<tReal>::operator() (tRand& aRng, unsigned aTid, GlobalData<tReal>&, tRandData& aEngData ) ->result_type
	{
		constexpr std::size_t digits_ = std::numeric_limits<tReal>::digits;
		return generate_canonical<tReal,digits_>( aRng, aTid, aEngData );
	}

	template< typename tReal > __host__ inline
	auto UniformReal::initialize( Identity<tReal>, std::size_t ) -> GlobalData<tReal>
	{
		return {};
	}

	template< typename tReal > __host__ inline
	void UniformReal::cleanup( GlobalData<tReal>& )
	{}
}
