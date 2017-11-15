#ifndef HISTOGRAM2D_CUH_8A47873B_D52D_43A0_8D3F_1E368AB2B386
#define HISTOGRAM2D_CUH_8A47873B_D52D_43A0_8D3F_1E368AB2B386

#include <cstddef>

namespace cusim
{
	template< typename tCount >
	class Histogram2D
	{
		public:
			using Count = tCount;

		public:
			__host__ Histogram2D( std::size_t, std::size_t, Count const* = nullptr );
			__host__  ~Histogram2D();

			__host__ Histogram2D( Histogram2D const& ) = delete;
			__host__ Histogram2D& operator= (Histogram2D const&) = delete;

			__host__ Histogram2D( Histogram2D&& ) noexcept;
			__host__ Histogram2D& operator= (Histogram2D&&) noexcept;

		public:
			__host__ void clear();
			__host__ void clear_async( cudaStream_t = 0 );

			__host__ Count* device_ptr();
			__host__ Count const* device_ptr() const;

			__host__ std::size_t n() const;
			__host__ std::size_t m() const;

			__host__ void free_cuda_resources();

		private:
			std::size_t mN, mM;
			Count* mDevBuffer;
	};


	template< typename tCount, typename tReal >
	class HistogramRecorder
	{
		public:
			using value_type = tCount;
		
		public:
			__host__ 
			HistogramRecorder( Histogram2D<tCount>&, tReal aDDE );

		public:
			__device__
			void record( tCount, tReal );

		private:
			tReal mDDE;
			tCount mDEBinCount; // and also stride in the histogram
			tCount* mHistoBuffer;
	};
}

#include "histogram2d.inl"
#endif // HISTOGRAM2D_CUH_8A47873B_D52D_43A0_8D3F_1E368AB2B386
