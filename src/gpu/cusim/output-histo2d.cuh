#ifndef OUTPUT_HISTO2D_CUH_D6E19EC7_E304_4FD8_9DF6_B4F1471CB1CF
#define OUTPUT_HISTO2D_CUH_D6E19EC7_E304_4FD8_9DF6_B4F1471CB1CF

#include "../support/utility.hpp"

#error "DEPRECATED - REMOVE"

namespace cusim
{
	template< typename tCount, typename tReal >
	class OutputHisto2D
	{
		public:
			using Real = tReal;
			using Count = tCount;
	
			struct Recorder
			{
				using value_type = tCount;
				
				__device__
				void record( tCount, tReal );

				tCount deBinCount; // and stride
				tReal dDE;
				tCount* buffer;
			};
	
		public:
			__host__
			OutputHisto2D( std::size_t, std::size_t, tReal, tCount const* = nullptr );

			__host__
			~OutputHisto2D();

		public:
			__host__
			void reset();

			__host__
			Recorder recorder() const;

			//TODO: distance to other???
			

			Count* HACK_buffer() const { return mDevBuffer; }

		private:
			std::size_t mKs, mDEs;
			tReal mDDE;
			Count* mDevBuffer;
	};
};

#include "output-histo2d.inl"
#endif // OUTPUT_HISTO2D_CUH_D6E19EC7_E304_4FD8_9DF6_B4F1471CB1CF
