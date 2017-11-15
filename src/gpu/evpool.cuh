/*-******************************************************* -- HEADER -{{{1- */
/*-	evpool.cuh : simple cudaEvent_t pool allocator
 */
/*-***************************************************************** -}}}1- */

#ifndef EVPOOL_CUH_DFCA34F1_5156_42DB_B824_3487424E16ED
#define EVPOOL_CUH_DFCA34F1_5156_42DB_B824_3487424E16ED

#include <vector>
#include <cstddef>

class EvPool final
{
	public:
		explicit EvPool( std::size_t aInitial = 256 );
		~EvPool() noexcept;

		EvPool( EvPool const& ) = delete;
		EvPool& operator= (EvPool const&) = delete;

		EvPool( EvPool&& ) noexcept;
		EvPool& operator= (EvPool&&) noexcept;

	public:
		cudaEvent_t alloc();
		void free( cudaEvent_t );

		void free_cuda_resources();

	public:
		std::vector<cudaEvent_t> mFreePool;
		std::size_t mAllocated;
};

#include "detail/evpool.inl"
#endif // EVPOOL_CUH_DFCA34F1_5156_42DB_B824_3487424E16ED

