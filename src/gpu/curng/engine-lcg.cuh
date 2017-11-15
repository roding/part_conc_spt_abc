/*-******************************************************* -- HEADER -{{{1- */
/*-	Random Engine: LCG
 *
 * Linear congurential based random engines.
 */
/*-***************************************************************** -}}}1- */

#ifndef ENGINE_LCG_CUH_416F1A1A_0A73_4F69_ABF8_9D54A321E985
#define ENGINE_LCG_CUH_416F1A1A_0A73_4F69_ABF8_9D54A321E985

#include <cstddef>
#include <cstdint>

#include "../support/utility.hpp"

namespace curng
{
	/** Random engine: Linear congurential
	 *
	 * Linear congurential random engine, with a 48 bit internal state and 32
	 * bits of output per invocation.
	 *
	 * TODO: reference; this is a POSIX setup, [jm]rand48
	 */
	struct EngineLCG48_32
	{
		struct GlobalData
		{
			using Factory = EngineLCG48_32;
			
			std::uint64_t* states;
		};

		class Engine
		{
			public:
				using result_type = std::uint32_t;

				static constexpr std::size_t bits = 32;

			public:
				__device__
				Engine( unsigned, GlobalData const& );

			private:
				static constexpr std::uint64_t a = 25214903917ull;
				static constexpr std::uint64_t c = 11ull;

			public:
				__device__
				result_type operator() (unsigned, GlobalData&);

			public:
				__device__
				void store( unsigned, GlobalData& );

				static __host__ __device__ constexpr
				result_type min(), max();

			private:
				std::uint64_t mState;
		};

		template< class tHostRNG > static __host__
		GlobalData initialize( std::size_t, tHostRNG& );

		static __host__
		void cleanup( GlobalData& );
	};
}

#include "engine-lcg.inl"
#endif // ENGINE_LCG_CUH_416F1A1A_0A73_4F69_ABF8_9D54A321E985
