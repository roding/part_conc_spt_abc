#ifndef DSINLINE_CUH_30B4BE51_206B_42F4_BC95_7A55E2CD0DFB
#define DSINLINE_CUH_30B4BE51_206B_42F4_BC95_7A55E2CD0DFB

#include <type_traits>
#include "../support/dsvalue.hpp"

namespace cusim
{
	/** Possibly inline storage
	 *
	 * `DSInline` facilitates passing data to a CUDA kernel via the kernel's
	 * function arguments. Depending on the element count \a tCount, it holds
	 * one of the three following three:
	 *   - a pointer to `tType` (dynamic element count)
	 *   - an array of tType (static element count)
	 *   - a single tType value (static element count of one)
	 *
	 * Static element counts require no further allocations; data is passed by
	 * value into the kernel (and copies are made by the CUDA runtime behind
	 * the scenes). Dynamic element counts require the user to allocate buffers
	 * in GPU memory (and to populate these).
	 *
	 * A number of helper functions are provided. See
	 *  - `acquire_host_ptr()`
	 *  - `release_host_ptr()`
	 *  - `upload_from_host_ptr()`
	 *  - `clean_gpu_cache()`
	 *
	 * Example:
	 * \code
	 * template< class tElementCount >
	 * struct MyKernelArguments
	 * {
	 * 	 ...
	 * 	 DSInline<float,tElementCount> elements;
	 *   ...
	 * };
	 *
	 * template< class tElementCount > 
	 * __global__ kernel( MyKernelArguments<tElementCount> aArguments )
	 * {
	 *   ...
	 * }
	 * \endcode
	 *
	 * See the methods `acquire_host_ptr()` and `upload_from_host_ptr()` for
	 * additional usage examples.
	 *
	 * \see acquire_host_ptr(), release_host_ptr(), upload_from_host_ptr(),
	 * clean_gpu_cache()
	 */
	template< typename tType, class tCount >
	struct DSInline
	{
		tType* values;

		static_assert( 
			std::is_same<typename tCount::type,unsigned>::value,
			"cusim: use unsigned counts"
		);
	};

	template< typename tType, unsigned tCount >
	struct DSInline<tType, StaticValue<unsigned,tCount>>
	{
		tType values[tCount];
	};

	template< typename tType >
	struct DSInline<tType, StaticValue<unsigned,1>>
	{
		tType value;
	};

	/** Acquire a host pointer to `DSInline` values
	 *
	 * `acquire_host_ptr()` returns a host-writable pointer, for the purpose
	 * of setting the values of a `DSInline` instance. The count \a tCount
	 * is either a `StaticValue<unsigned,VALUE>` or a `DynamicValue<unsigned>`.
	   \code
	   using Count_ = ...; // StaticValue or DynamicValue
	   Count_ count{ ... };
	  
	   ...
	   DSInline<float,Count_> dsi;
	   float* hptr = acquire_host_ptr( dsi, count );

	   ...
	   while( running ) {
	     hptr[0] = 44.f;
	     ...
	  
	     upload_from_host_ptr( dsi, hptr, count, ... );
	     kernel<<<...>>>( dsi );
	     clean_gpu_cache( dsi, ... );
	   }
	  
	   ...
	   release_host_ptr( dsi, hptr );
	   \endcode
	 *
	 * Pointers returned by `acquire_host_ptr()` must later be released with
	 * `release_host_ptr()`.
	 *
	 * <i>Internal details:</i> For static element counts, `acquire_host_ptr()`
	 * simply returns a pointer to the internal (inlined) storage of the 
	 * `DSInline` object. For dynamic element counts, `acquire_host_ptr()` will
	 * allocate (via `new []`) host storage.
	 */
	template< typename tType, class tCount > __host__
	tType* acquire_host_ptr( DSInline<tType,tCount>&, tCount const& );

	/** Release host pointer of `DSInline`
	 *
	 * Release a host pointer previously acquired using `acquire_host_ptr()`.
	 * The host pointer will no longer be usable.
	 *
	 * <i>Internal details:</i> For static element counts, `release_host_ptr()`
	 * does nothing. For dynamic element counts, `release_host_ptr()` 
	 * deallocates (via `delete []`) the pointer.
	 */
	template< typename tType, class tCount > __host__
	void release_host_ptr( DSInline<tType,tCount>&, tType* );

	/** Upload data from host storage to `DSInline` object
	 *
	 * Uploads \a tCount values from host storage to the `DSInline` object.
	 *
	 * \warning The host storage pointer `tType const*` <b>must</b> have been
	 * returned by `acquire_host_ptr()`. It is not permissible to use any other
	 * host pointers.
	 *
	 * If necessary, `upload_from_host_ptr()` will attempt to acquire GPU
	 * storage via the \a tAlloc callback. When called, \a tAlloc must return
	 * a suitable GPU memory pointer (i.e., it must be sized sufficiently to
	 * store \a tCount values of \a tType).
	 *
	 * Example:
	   \code
	   Count_ count{ ... };
	   DSInline<float,Count_> dsi;

	   ...
	   float* hptr = acquire_host_ptr( dsi, count );

	   ...
	   upload_from_host_ptr(
	   	dsi,
	   	hptr,
	   	count,
	   	[&] () { float* ret; cudaMalloc( &ret, sizeof(float)*count ); return ret; }
	   );
	  
	   ...
	   clean_gpu_cache( dsi, [&] (float* aPtr) { cudaFree( aPtr ); } );
	   
	   ...
	   release_host_ptr( dsi, hptr );
	   \endcode
	 *
	 * Each `upload_from_host_ptr()` call must be matched with a call to
	 * `clean_gpu_cache()`.
	 * 
	 * <i>Internal details:</i> For static element counts,
	 * `upload_from_host_ptr()` does nothing (\a tAlloc is not invoked!). For
	 * dynamic element counts, `upload_from_host_ptr()` acquires GPU storage
	 * via \a tAlloc, stores the pointer to the storage in the \a DSInline
	 * object, and finally issues a `cudaMemcpyAsync()` from the host storage
	 * to the acquired GPU storage.
	 */
	template< typename tType, class tCount, class tAlloc > __host__
	void upload_from_host_ptr( 
		DSInline<tType,tCount>&, 
		tType const*, 
		tCount const&, 
		tAlloc&&,
		cudaStream_t = 0
	);

	/** Clean temporary GPU objects
	 *
	 * Clean up temporary GPU objects from a call to `upload_from_host_ptr()`.
	 * Data stored in the \a DSInline object may become unavailable; care must
	 * be taked to only call `clean_gpu_cache()` after all kernels using the
	 * `DSInline` object have finished.
	 *
	 * <i>Internal details:</i> For static element counts, `clean_gpu_cache()`
	 * does nothing (and \a tFree is never invoked!). For dynamic element
	 * counts, `clean_gpu_cache()` releases the GPU storage via a call to \a
	 * tFree.
	 */
	template< typename tType, class tCount, class tFree > __host__
	void clean_gpu_cache( DSInline<tType,tCount>&, tFree&& );
}

#include "dsinline.inl"
#endif // DSINLINE_CUH_30B4BE51_206B_42F4_BC95_7A55E2CD0DFB
