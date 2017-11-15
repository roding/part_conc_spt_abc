#include "common.hpp"

#include <sstream>
#include <stdexcept>

#include <cstring>

#include <tinyformat/tinyformat.h>

#include "errors.hpp"

std::string to_string( EModel aModel )
{
	switch( aModel )
	{
		case EModel::discreteFixedZ: return "discrete-fixed-depth";
		case EModel::discreteVariableZ: return "discrete-variable-depth";
		
		case EModel::lognormal: return "lognormal";
	}

	return tfm::format( "EModel(%ld)", long(aModel) );
}
std::string to_string( EWeightingScheme aScheme )
{
	switch( aScheme )
	{
		case EWeightingScheme::pmcStandard: return "pmc-standard";
		case EWeightingScheme::inverseDistSq: return "inverse-distance-squared";
	}


	return tfm::format( "EWeightingScheme(%ld)", long(aScheme) );
}

template<>
EModel from_string( char const* aStr, Identity<EModel> )
{
	if( 0 == std::strcmp( "discrete-fixed-depth", aStr ) )
		return EModel::discreteFixedZ;
	else if( 0 == std::strcmp( "discrete-variable-depth", aStr ) )
		return EModel::discreteVariableZ;
	else if( 0 == std::strcmp( "discrete", aStr ) )
	{
		std::fprintf( stderr, "WARNING: EModel: deprecated value \"discrete\"\n" );
		return EModel::discrete;
	}
	else if( 0 == std::strcmp( "lognormal", aStr ) )
	{
		std::fprintf( stderr, "WARNING: EModel: deprecated value \"lognormal\"\n" );
		return EModel::lognormal;
	}

	long dist = 0;
	if( 1 == std::sscanf( aStr, "EModel(%ld)", &dist ) )
		return EModel(dist);

	throw error::EnumUnknownString( "'%s' is not a valid EModel", aStr );
}
template<>
EWeightingScheme from_string( char const* aStr, Identity<EWeightingScheme> )
{
	if( 0 == std::strcmp( "pmc-standard", aStr ) )
		return EWeightingScheme::pmcStandard;
	else if( 0 == std::strcmp( "inverse-distance-squared", aStr ) )
		return EWeightingScheme::inverseDistSq;
	

	long scheme = 0;
	if( 1 == std::sscanf( aStr, "EWeightingScheme(%ld)", &scheme ) )
		return EWeightingScheme(scheme);

	throw error::EnumUnknownString( "'%s' is not a valid EWeightingScheme", aStr );
}
