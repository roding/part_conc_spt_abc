#include <limits>
#include <type_traits>

#include "support.cuh"

namespace curng
{
	namespace detail
	{
		__device__ inline
		void sincos_(float aX, float* aSin, float* aCos)
		{
			sincosf(aX, aSin, aCos);
		}
		__device__ inline
		void sincos_(double aX, double* aSin, double* aCos)
		{
			sincos(aX, aSin, aCos);
		}
	}


	template< typename tReal > __device__ inline
	NormalBoxMuller::Distribution<tReal>::Distribution( unsigned, NormalBoxMuller::GlobalData<tReal> const& )
		: mCache( std::numeric_limits<tReal>::quiet_NaN() )
	{}

	template< typename tReal > template< class tRand, class tRandData > __device__ inline
	auto NormalBoxMuller::Distribution<tReal>::operator() (tRand& aRng, unsigned aTid, GlobalData<tReal>&, tRandData& aEngData ) -> result_type
	{
		constexpr tReal pi2 = tReal(2)*tReal(3.14159265357989);
		constexpr tReal eps = std::numeric_limits<tReal>::min();
		constexpr std::size_t digits_ = std::numeric_limits<tReal>::digits;

		if( !/*std::*/isnan( mCache ) )
		{
			auto ret = mCache;
			mCache = std::numeric_limits<tReal>::quiet_NaN();
			return ret;
		}

		tReal u1;
		do
		{
			u1 = generate_canonical<tReal,digits_>( aRng, aTid, aEngData );
		} while (u1 <= eps);

		tReal const lu2 = std::sqrt(tReal(-2) * std::log(u1));
		
		tReal u2 = generate_canonical<tReal,digits_>( aRng, aTid, aEngData );

		tReal s, c;
		detail::sincos_(pi2 * u2, &s, &c);

		mCache = lu2 * s;
		return lu2 * c;
	}

	template< typename tReal > __host__ inline
	auto NormalBoxMuller::initialize( Identity<tReal>, std::size_t ) -> GlobalData<tReal>
	{
		return {};
	}

	template< typename tReal > __host__ inline
	void NormalBoxMuller::cleanup( GlobalData<tReal>& )
	{}
}
