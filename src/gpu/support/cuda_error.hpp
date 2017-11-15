#ifndef CUDA_ERROR_HPP_AF1800D2_0D52_4782_B694_CF9BF1FC8CA7
#define CUDA_ERROR_HPP_AF1800D2_0D52_4782_B694_CF9BF1FC8CA7

#include <system_error>

#include <cuda_runtime.h>

#define CUDA_CHECKED ::detail::CudaErrorChecker(__FILE__,__LINE__,__func__) = 


namespace error
{
	struct CudaError
		: std::system_error
	{
		using std::system_error::system_error;
	};
}


std::error_category const& cuda_error_category() noexcept;

std::error_code make_cuda_error_code() noexcept;
std::error_code make_cuda_error_code( cudaError_t ) noexcept;


namespace detail
{
	struct CudaErrorChecker
	{
		CudaErrorChecker( char const*, int, char const* );

		cudaError_t operator= (cudaError_t) const;

		char const* file;
		char const* func;
		int line;
	};
}

#endif // CUDA_ERROR_HPP_AF1800D2_0D52_4782_B694_CF9BF1FC8CA7


