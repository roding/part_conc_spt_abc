#include "utility.hpp"

#include <stdexcept>

namespace detail
{
	[[noreturn]] void throw_bad_int_cast()
	{
		throw std::logic_error( "Bad integer cast" );
	}
}
