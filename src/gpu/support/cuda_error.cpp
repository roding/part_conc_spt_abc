#include "cuda_error.hpp"

#include <tinyformat/tinyformat.h>

namespace
{
	struct CudaErrorCategory_ : std::error_category
	{
		char const* name() const noexcept override;
		std::string message( int ) const override;
	};
}


std::error_category const& cuda_error_category() noexcept
{
	static CudaErrorCategory_ cudaCat;
	return cudaCat;
}

std::error_code make_cuda_error_code() noexcept
{
	return std::error_code( cudaSuccess, cuda_error_category() );
}
std::error_code make_cuda_error_code( cudaError_t aErrCode ) noexcept
{
	return std::error_code( aErrCode, cuda_error_category() );
}


namespace detail
{
	CudaErrorChecker::CudaErrorChecker( char const* aFile, int aLine, char const* aFunc )
		: file(aFile)
		, func(aFunc)
		, line(aLine)
	{}

	cudaError_t CudaErrorChecker::operator= (cudaError_t aErrCode) const
	{
		if( cudaSuccess != aErrCode )
		{
			throw error::CudaError(
				make_cuda_error_code( aErrCode ),
				tfm::format( "from %s() at %s:%d", func, file, line )
			);
		}

		return aErrCode;
	}
}

namespace
{
	char const* CudaErrorCategory_::name() const noexcept
	{
		return "CUDA Error";
	}

	std::string CudaErrorCategory_::message( int aErrCode ) const
	{
		cudaError_t const ec = cudaError_t(aErrCode);
		return tfm::format( "%s (%s)", cudaGetErrorString(ec), cudaGetErrorName(ec) );
	}
}
