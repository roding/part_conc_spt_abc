#include <limits>
#include <cstdint>

namespace curng
{
	namespace detail
	{
		__host__ __device__ constexpr
		unsigned log2_( std::uint64_t aX )
		{
			/* The special case of aX == 0 is intentional (and emulates the
			 * behaviour of the __log2<> in libc++. Essentially, it's required
			 * for the case where we use log2_() on UINT64_MAX-UIN64_MIN+1
			 * (which is equal to zero due to overflow) below. In this case,
			 * log2(0) will return 64.
			 */
			return (aX==0) 
				? 64
				: (aX<=1) ? 0 : 1+log2_(aX/2)
			;
		}
	}

	
	template< typename tReal, std::size_t tBits, class tEngine, class tEngData >
	__device__ inline
	tReal generate_canonical( tEngine& aEng, unsigned aTid, tEngData& aEngData )
	{
		constexpr unsigned kDigits_ = std::numeric_limits<tReal>::digits;

		constexpr unsigned need_ = kDigits_ < tBits ? kDigits_ : tBits;
		constexpr unsigned get_ = detail::log2_(tEngine::max()-tEngine::min()+1);
		constexpr unsigned rounds_ = need_/get_ + (need_%get_ != 0) + (need_==0);

		constexpr tReal mul_ = tEngine::max() - tEngine::min() + tReal(1);

#		if 0
		/* Broken version, returns [0,1] inclusive */
		tReal base = mul_;
		tReal acc = aEng( aTid, aEngData ) - tEngine::min();
		for( unsigned i = 1; i < rounds_; ++i, base *= mul_ )
			acc += (aEng( aTid, aEngData ) - tEngine::min()) * base;

		return acc/base;
#		else
		/* "Fixed" version, rejects numbers >= 1 */
		tReal acc, base;
		do
		{
			base = mul_;
			acc = aEng( aTid, aEngData ) - tEngine::min();
			for( unsigned i = 1; i < rounds_; ++i, base *= mul_ )
				acc += (aEng( aTid, aEngData ) - tEngine::min()) * base;
		} while( acc >= base );

		return acc/base;
#		endif
	}
}
