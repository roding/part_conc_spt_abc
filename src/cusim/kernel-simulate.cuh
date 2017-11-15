#ifndef KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009
#define KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009

namespace cusim
{
	template< 
		class tSystemSetup,
		class tRunData,
		class tOutput,
		class tRandData, 
		class tReal = typename tSystemSetup::value_type,
		class tCount = typename tSystemSetup::count_type,
		class tRandomFactory = typename tRandData::Factory
	>
	__global__
	void K_simulate_system( tSystemSetup, tRunData, tOutput, tRandData );
}

#include "kernel-simulate.inl"
#endif // KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009
