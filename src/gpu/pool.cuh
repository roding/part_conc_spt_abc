/*-******************************************************* -- HEADER -{{{1- */
/*-	pool.cuh : simple device memory pool
 */
/*-***************************************************************** -}}}1- */

#ifndef POOL_CUH_3E195F29_A244_40EE_AF6C_5F3B3F0AE82F
#define POOL_CUH_3E195F29_A244_40EE_AF6C_5F3B3F0AE82F

#include <tuple>
#include <vector>

#include <cstddef>

/** Device memory pool for array allocations of \a tType
 *
 * `Pool<>` provides any number of array allocations of \a tType; the size of
 * each individual array is selected at runtime (but is fixed for the lifetime
 * of a `Pool<>` instance). \a tType should probably be trivially constructible
 * and destructible. 
 *
 * `Pool<>` is useful when running multiple simulations asynchronously in
 * parallel.
 */
template< typename tType >
class Pool final
{
	public:
		/** Constructor: `Pool<>`
		 *
		 * Constructs a `Pool<>` where each returned device memory allocation
		 * is an array `tType[aArrayCount]`. The \a aInitial parameter controls
		 * the number of initial allocations that the pool holds (the pool can
		 * grow at runtime, though).
		 */
		explicit Pool( std::size_t aArrayCount, std::size_t aInitial = 128 );
		~Pool() noexcept;

		Pool( Pool const& ) = delete;
		Pool& operator= (Pool const&) = delete;

		Pool( Pool&& ) noexcept;
		Pool& operator= (Pool&&) noexcept;

	public:
		/** Return single allocation
		 *
		 * Returns a device pointer to an allocation of `tType[arrayCount]`,
		 * where \a arrayCount was specified during construction. Each 
		 * allocation from `alloc()` must be freed with a matching `free()`.
		 */
		tType* alloc();
		/** Free allocation
		 * 
		 * Free an allocation previously allocated with `alloc()`.
		 */
		void free( tType const* );

		/** Free CUDA resources
		 *
		 * Free CUDA resources held by this pool. This ensures that no CUDA 
		 * API functions need to be called by the destructor.
		 */
		void free_cuda_resources();

	public:
		std::vector<tType*> mFreePool;
		std::size_t mArrayCount;
		std::size_t mAllocated;
};

/** Device-mapped host memory pool for array allocations of \a tType
 *
 * See `Pool<>` for general information. `MappedPool` allocates device-mappable
 * host memory instead of dedicated device memory. 
 *
 * This allows the distance computation to write its result directly into
 * host memory, removing the need for a 4-byte device-to-host asynchronous
 * memcopy operation.
 */
template< typename tType >
class MappedPool final
{
	public:
		/** Constructor: `Pool<>`
		 *
		 * Constructs a `Pool<>` where each returned device memory allocation
		 * is an array `tType[aArrayCount]`. The \a aInitial parameter controls
		 * the number of initial allocations that the pool holds (the pool can
		 * grow at runtime, though).
		 */
		explicit MappedPool( std::size_t aArrayCount, std::size_t aInitial = 128 );
		~MappedPool() noexcept;

		MappedPool( MappedPool const& ) = delete;
		MappedPool& operator= (MappedPool const&) = delete;

		MappedPool( MappedPool&& ) noexcept;
		MappedPool& operator= (MappedPool&&) noexcept;

	public:
		/** Return single allocation
		 *
		 * Returns a host and device pointer, both refering to the same 
		 * memory.
		 */
		std::tuple<tType*,tType*> alloc();
		/** Free allocation
		 * 
		 * Free an allocation previously allocated with `alloc()`.
		 */
		void free( tType const* aHost, tType const* aDev );

		/** Free CUDA resources
		 *
		 * Free CUDA resources held by this pool. This ensures that no CUDA 
		 * API functions need to be called by the destructor.
		 */
		void free_cuda_resources();

	public:
		std::vector<tType*> mFreePoolHost;
		std::vector<tType*> mFreePoolDev;
		std::size_t mArrayCount;
		std::size_t mAllocated;
};


#include "detail/pool.inl"
#endif // POOL_CUH_3E195F29_A244_40EE_AF6C_5F3B3F0AE82F

