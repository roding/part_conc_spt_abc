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
	
	template< class... tFixed >
	SimPtr_ finalize_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		return make_unique<SimulationT<
			tFixed...,
			sim_arg::HostRng<SimHostRNG>
		>>( aPar, aRng, aCfg );
	}

	template< class... tFixed >
	SimPtr_ fix_distance_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.distanceType )
		{
			case ESimScalar::floatType: 
				return finalize_< tFixed..., sim_arg::DistanceType<float> >(aRng,aPar,aCfg);
			case ESimScalar::doubleType:
				return finalize_< tFixed..., sim_arg::DistanceType<double> >(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Distance type: scalar type ESimScalar(%ld) not supported", long(aCfg.deviceScalarType) ) );
	}
	template< class... tFixed >
	SimPtr_ fix_count_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.countType )
		{
			case ESimCount::uint32: 
				return fix_distance_<tFixed..., sim_arg::CountType<std::uint32_t>>(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Count type ESimCount(%ld) not handled", long(aCfg.countType) ) );
	}

	template< class... tFixed >
	SimPtr_ fix_dscalar_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.deviceScalarType )
		{
			case ESimScalar::floatType: 
				return fix_count_< tFixed..., sim_arg::DeviceScalar<float> >(aRng,aPar,aCfg);
			case ESimScalar::doubleType:
				return fix_count_< tFixed..., sim_arg::DeviceScalar<double> >(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Device scalar type: scalar type ESimScalar(%ld) not supported", long(aCfg.deviceScalarType) ) );
	}
	template< class... tFixed >
	SimPtr_ fix_hscalar_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aCfg.hostScalarType )
		{
			case ESimScalar::floatType: 
				return fix_dscalar_< tFixed..., sim_arg::HostScalar<float> >(aRng,aPar,aCfg);
			case ESimScalar::doubleType:
				return fix_dscalar_< tFixed..., sim_arg::HostScalar<double> >(aRng,aPar,aCfg);
		}

		throw std::runtime_error( tfm::format( "Host scalar type: scalar type ESimScalar(%ld) not supported", long(aCfg.hostScalarType) ) );
	}

	template< class... tFixed >
	SimPtr_ fix_ccount_( SimHostRNG& aRng, Param_ aPar, Conf_ aCfg )
	{
		switch( aPar.componentCount )
		{
			case 1: return fix_hscalar_< tFixed..., sim_arg::ComponentCount<SVal_<1>> >( aRng, aPar, aCfg );
			case 2: return fix_hscalar_< tFixed..., sim_arg::ComponentCount<SVal_<2>> >( aRng, aPar, aCfg );
		}
		
		return fix_hscalar_<tFixed..., sim_arg::ComponentCount<DVal_> >( aRng, aPar, aCfg );
	}
}

std::unique_ptr<Simulation> make_simulation( SimHostRNG& aRng, input::Parameters const& aPar, SimulationConfig const& aCfg )
{
	// Fix the EModel
	switch( aPar.model )
	{
		case EModel::discreteFixedZ: 
			return fix_ccount_< sim_arg::Model<EModel::discreteFixedZ> >( aRng, aPar, aCfg );

		case EModel::discreteVariableZ: 
			return fix_ccount_< sim_arg::Model<EModel::discreteVariableZ> >( aRng, aPar, aCfg );
	}

	throw std::runtime_error( tfm::format( "make_simulation(): EModel %s not handled\n", to_string(aPar.model) ) );
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
