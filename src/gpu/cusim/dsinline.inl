#include "../support/cuda_error.hpp"

namespace cusim
{
	namespace detail
	{
		template< class tCount >
		struct DSOps_
		{
			template< typename tType > static inline
			tType* acquire( DSInline<tType,tCount>&, tCount const& aCount )
			{
				return new tType[aCount];
			}
			template< typename tType > static inline
			void release( DSInline<tType,tCount>&, tType* aHostPtr )
			{
				delete [] aHostPtr;
			}

			template< typename tType, class tAlloc > static inline
			void upload( DSInline<tType,tCount>& aInline, tType const* aHostPtr, tCount const& aCount, cudaStream_t aStream, tAlloc&& aAlloc )
			{
				aInline.values = std::forward<tAlloc>(aAlloc)();

				CUDA_CHECKED cudaMemcpyAsync(
					aInline.values,
					aHostPtr,
					sizeof(tType)*aCount,
					cudaMemcpyHostToDevice,
					aStream
				);
			}
			template< typename tType, class tFree > static inline
			void clean( DSInline<tType,tCount>& aInline, tFree&& aFree )
			{
				std::forward<tFree>(aFree)( aInline.values );
			}
		};

		template< unsigned tCount >
		struct DSOps_< StaticValue<unsigned,tCount> >
		{
			using Count_ = StaticValue<unsigned,tCount>;
			
			template< typename tType > static inline
			tType* acquire( DSInline<tType,Count_>& aInline, Count_ const& )
			{
				return aInline.values;
			}
			template< typename tType > static inline
			void release( DSInline<tType,Count_>&, tType* )
			{}

			template< typename tType, class tAlloc > static inline
			void upload( DSInline<tType,Count_>&, tType const*, Count_ const&, cudaStream_t, tAlloc&& )
			{}
			template< typename tType, class tFree > static inline
			void clean( DSInline<tType,Count_>&, tFree&&  )
			{}
		};
		template<>
		struct DSOps_< StaticValue<unsigned,1> >
		{
			using Count_ = StaticValue<unsigned,1>;

			template< typename tType > static inline
			tType* acquire( DSInline<tType,Count_>& aInline, Count_ const& )
			{
				return &aInline.value;
			}
			template< typename tType > static inline
			void release( DSInline<tType,Count_>&, tType* )
			{}

			template< typename tType, class tAlloc > static inline
			void upload( DSInline<tType,Count_>&, tType const*, Count_ const&, cudaStream_t, tAlloc&& )
			{}
			template< typename tType, class tFree > static inline
			void clean( DSInline<tType,Count_>&, tFree&&  )
			{}
		};

	}

	template< typename tType, class tCount > __host__ inline
	tType* acquire_host_ptr( DSInline<tType,tCount>& aInline, tCount const& aCount )
	{
		return detail::DSOps_<tCount>::acquire( aInline, aCount );
	}

	template< typename tType, class tCount > __host__ inline
	void release_host_ptr( DSInline<tType,tCount>& aInline, tType* aHostPtr )
	{
		detail::DSOps_<tCount>::release( aInline, aHostPtr );
	}

	template< typename tType, class tCount, class tAlloc > __host__ inline
	void upload_from_host_ptr( DSInline<tType,tCount>& aInline, tType const* aHostPtr, tCount const& aCount, tAlloc&& aAlloc, cudaStream_t aStream )
	{
		detail::DSOps_<tCount>::upload( aInline, aHostPtr, aCount, aStream, std::forward<tAlloc>(aAlloc) );
	}

	template< typename tType, class tCount, class tFree > __host__ inline
	void clean_gpu_cache( DSInline<tType,tCount>& aInline, tFree&& aFree )
	{
		detail::DSOps_<tCount>::clean( aInline, std::forward<tFree>(aFree) );
	}
}
