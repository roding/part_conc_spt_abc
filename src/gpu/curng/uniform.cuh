/*-******************************************************* -- HEADER -{{{1- */
/*-	Uniform distribution
 */
/*-***************************************************************** -}}}1- */

#ifndef UNIFORM_CUH_C74AA5EB_0AC7_45B6_8A34_4CD8D76CC927
#define UNIFORM_CUH_C74AA5EB_0AC7_45B6_8A34_4CD8D76CC927

#include "support.cuh"

#include "../support/utility.hpp"

namespace curng
{
	/** Uniform real distribution
	 *
	 * Produce uniformly distributed real numbers. Uses an interface similar
	 * to the normal distributions (e.g., NormalBoxMuller). Essentially just
	 * returns `curng::generate_canonical<...>(), but with a more complicated
	 * interface.
	 */
	struct UniformReal 
	{
		template< typename tReal >
		struct GlobalData
		{
			using Factory = UniformReal;
			using result_type = tReal;
		};

		template< typename tReal >
		class Distribution
		{
			public:
				using result_type = tReal;

			public:
				__device__
				Distribution( unsigned, GlobalData<tReal> const& );

			public:
				template< class tRand, class tRandData > __device__
				result_type operator() (tRand&, unsigned, GlobalData<tReal>&, tRandData& );
		};

		template< typename tReal > static __host__
		GlobalData<tReal> initialize( Identity<tReal>, std::size_t );

		template< typename tReal > static __host__
		void cleanup( GlobalData<tReal>& );
	};
}

#include "uniform.inl"
#endif // UNIFORM_CUH_C74AA5EB_0AC7_45B6_8A34_4CD8D76CC927
