#include "simulation_impl.cuh"

#include <stdexcept>

#include <cstddef>
#include <cstdint>

#include <tinyformat/tinyformat.h>

#include "shared/errors.hpp"
#include "support/utility.hpp"

#include "simulation.hpp"

/* TODO: make a few specific instantiations compile in different .cpp
 *
 * Note: for extern foo<bar,baz> we'll probably have to match the order of the 
 * template declaration. An alternative is to make make finalize_() external.
 */

namespace
{
	using Param_ = input::Parameters const&;
	using Conf_ = SimulationConfig const&;
	using SimPtr_ = std::unique_ptr<Simulation>;

	using DVal_ = DynamicValue<unsigned>;
	template< unsigned tVal >
	using SVal_ = StaticValue<unsigned,tVal>;
	
	template< EModel tModel, class tCCount, typename tScalar, typename tCount >
	SimPtr_ finalize_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		return make_unique<SimulationT<
			sim_arg::ScalarType<tScalar>,
			sim_arg::CountType<tCount>,
			sim_arg::ComponentCount<tCCount>,
			sim_arg::Model<tModel>,
			sim_arg::HostRng<SimHostRNG>
		>>( aPar, aRng, aCfg );
	}

	template< EModel tModel, class tCCount, typename tScalar >
	SimPtr_ fix_ctype_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.countType )
		{
			case ESimCount::uint32: 
				return finalize_<tModel,tCCount,float,std::uint32_t>(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Count type ESimCount(%ld) not handled", long(aCfg.countType) ) );
	}

	template< EModel tModel, class tCCount >
	SimPtr_ fix_count_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.scalarType )
		{
			case ESimScalar::floatType: 
				return fix_ctype_<tModel,tCCount,float>(aRng,aPar,aCfg);
			case ESimScalar::doubleType:
				return fix_ctype_<tModel,tCCount,double>(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Scalar type ESimScalar(%ld) not handled", long(aCfg.scalarType) ) );
	}
	template< EModel tModel >
	SimPtr_ fix_model_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aPar.componentCount )
		{
			case 1: return fix_count_<tModel, SVal_<1>>( aRng, aPar, aCfg );
			case 2: return fix_count_<tModel, SVal_<2>>( aRng, aPar, aCfg );
		}
		
		return fix_count_<tModel, DVal_>( aRng, aPar, aCfg );
	}
}

std::unique_ptr<Simulation> make_simulation( SimHostRNG& aRng, input::Parameters const& aPar, SimulationConfig const& aCfg )
{
	// Fix the distribution class
	switch( aPar.model )
	{
		case EModel::discreteFixedZ: 
			return fix_model_<EModel::discreteFixedZ>( aRng, aPar, aCfg );

		case EModel::discreteVariableZ: 
			return fix_model_<EModel::discreteVariableZ>( aRng, aPar, aCfg );
	}

	throw std::runtime_error( tfm::format( "make_simulation(): Model %s not handled\n", to_string(aPar.model) ) );
}


namespace detail
{
	[[noreturn]] void throw_logic_error_( char const* aErrMsg )
	{
		throw std::logic_error( aErrMsg ); //XXX-named error
	}
	[[noreturn]] void throw_invalid_gpuspec_( char const* aErrMsg )
	{
		throw error::InvalidGPUSpec( aErrMsg );
	}
}
