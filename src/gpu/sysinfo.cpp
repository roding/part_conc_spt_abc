#include "sysinfo.hpp"

#include <cstdio>
#include <cuda_runtime.h>

#include "support/cuda_error.hpp"

void list_devices()
{
	int count = 0;
	CUDA_CHECKED cudaGetDeviceCount( &count );
	
	std::printf( "%d devices:\n", count );
	for( int dev = 0; dev < count; ++dev )
	{
		cudaDeviceProp dp;
		CUDA_CHECKED cudaGetDeviceProperties( &dp, dev );

		char const* computeMode = [&] {
			switch( dp.computeMode )
			{
				case cudaComputeModeDefault: return "default";
				case cudaComputeModeExclusive: return "exclusive";
				case cudaComputeModeProhibited: return "prohibited";
				case cudaComputeModeExclusiveProcess: return "exclusive-process";
			}

			return "<unknown>";
		}();

		std::printf( " ID = %d.   \"%s\" (Compute %d.%d; %s; timeout %s)\n", dev+1, dp.name, dp.major, dp.minor, computeMode, dp.kernelExecTimeoutEnabled ? "yes" : "no" );

		std::printf( "   + SMs: %u (max conc. thread count = %u)\n", dp.multiProcessorCount, dp.multiProcessorCount*dp.maxThreadsPerMultiProcessor );
		std::printf( "   + warps per block: %u, per SM: %u\n", dp.maxThreadsPerBlock/32, dp.maxThreadsPerMultiProcessor/32 );
		std::printf( "   + smem per block: %zu kBytes, per SM: %zu kBytes\n", dp.sharedMemPerBlock/1024, dp.sharedMemPerMultiprocessor/1024 );
		std::printf( "   + concurrent kernels: %s, concurrent copy: %s, async eng: %u\n", dp.concurrentKernels ? "yes" : "no", dp.deviceOverlap ? "yes" : "no", dp.asyncEngineCount );
		std::printf( "   + can map hostmem: %s, unified addr: %s, pageable mem: %s\n", dp.canMapHostMemory ? "yes" : "no", dp.unifiedAddressing ? "yes" : "no", dp.pageableMemoryAccess ? "yes" : "no" );
		std::printf( "   + float to double perf. ratio: %dx\n", dp.singleToDoublePrecisionPerfRatio );
	}
}
