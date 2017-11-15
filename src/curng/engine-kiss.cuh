/*-******************************************************* -- HEADER -{{{1- */
/*-	Random Engine: KISS
 */
/*-***************************************************************** -}}}1- */

#ifndef ENGINE_KISS_CUH_C6E09780_42EF_4D51_A8F0_CEB8CCE2D05C
#define ENGINE_KISS_CUH_C6E09780_42EF_4D51_A8F0_CEB8CCE2D05C

#include <cstddef>
#include <cstdint>

namespace curng
{
	/* Random engine: KISS
	 * 
	 * Implementation of the "KISS" random engine, as defined in the paper
	 * 
	 *  - "A Comment on the Implementation of the Ziggurat Method",
     *     Leong, Zhang, Lee, Luk and Villasenor; 2005
	 *
	 * See definition of the `KISS` macro in the paper. Produces 32 bits out
	 * output each invocation.
	 */
	struct EngineKISS
	{
		struct GlobalData
		{
			using Factory = EngineKISS;
			
			std::uint32_t* jsr;
			std::uint32_t* z;
			std::uint32_t* w;
			std::uint32_t* jcong;
		};


		class Engine
		{
			public:
				using result_type = std::uint32_t;
		
				static constexpr std::size_t bits = 32;
			
			public:
				__device__
				Engine( unsigned, GlobalData const& );

			public:
				__device__
				result_type operator() (unsigned, GlobalData&);

			public:
				__device__
				void store( unsigned, GlobalData& );

				static __host__ __device__ constexpr
				result_type min(), max();

			private:
				std::uint32_t jsr;
				std::uint32_t z, w;
				std::uint32_t jcong;
		};

		template< class tHostRNG > static __host__
		GlobalData initialize( std::size_t, tHostRNG& );

		static __host__
		void cleanup( GlobalData& );
	};
}

#include "engine-kiss.inl"
#endif // ENGINE_KISS_CUH_C6E09780_42EF_4D51_A8F0_CEB8CCE2D05C
