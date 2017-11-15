#ifndef COMMON_HPP_FDE5CB6E_BC0F_492D_970E_F571571F326A
#define COMMON_HPP_FDE5CB6E_BC0F_492D_970E_F571571F326A

#include <string>

#include "../support/compat.hpp"
#include "../support/utility.hpp"

enum class EModel
{
	discreteFixedZ,
	discreteVariableZ,

	lognormal ATTR_DEPRECATED,
	discrete ATTR_DEPRECATED = discreteFixedZ
};

using EDistributionClass ATTR_DEPRECATED = EModel;

enum class EWeightingScheme
{
	pmcStandard,
	inverseDistSq
};



std::string to_string( EModel );
std::string to_string( EWeightingScheme );

template< typename tEnum >
tEnum from_string( char const*, Identity<tEnum> = Identity<tEnum>{} );

#endif // COMMON_HPP_FDE5CB6E_BC0F_492D_970E_F571571F326A
