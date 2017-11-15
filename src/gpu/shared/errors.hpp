#ifndef ERRORS_HPP_79A39D1D_A7CF_4CC0_A118_1D243C3076A6
#define ERRORS_HPP_79A39D1D_A7CF_4CC0_A118_1D243C3076A6

#include <stdexcept>

namespace error
{
	struct RuntimeError : std::runtime_error
	{
		template< typename... tArgs >
		RuntimeError( char const*, tArgs&&... );
	};

		struct XMLLoadError : RuntimeError
		{
			using RuntimeError::RuntimeError;
		};

			struct XMLMissing : XMLLoadError
			{
				using XMLLoadError::XMLLoadError;
			};
			struct XMLParsing : XMLLoadError
			{
				using XMLLoadError::XMLLoadError;
			};
			struct XMLGeneric : XMLLoadError
			{
				using XMLLoadError::XMLLoadError;
			};

		struct EnumUnknownString : RuntimeError
		{
			using RuntimeError::RuntimeError;
		};

		struct InvalidGPUSpec : RuntimeError
		{
			using RuntimeError::RuntimeError;
		};
}

#include "errors.inl"
#endif // ERRORS_HPP_79A39D1D_A7CF_4CC0_A118_1D243C3076A6
