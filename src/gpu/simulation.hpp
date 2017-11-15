/*-******************************************************* -- HEADER -{{{1- */
/*-	simulation.hpp
 *
 * Interface to a simulation instance. See `make_simulation()`.
 */
/*-***************************************************************** -}}}1- */

#ifndef SIMULATION_HPP_37C3055A_87E4_48BB_B6B8_1A4D16482DFC
#define SIMULATION_HPP_37C3055A_87E4_48BB_B6B8_1A4D16482DFC

#include <memory>
#include <random>
#include <string>

#include "shared/input.hpp"
#include "shared/common.hpp"

using SimHostRNG = std::mt19937;

enum class ESimScalar
{
	floatType,
	doubleType
};
enum class ESimCount
{
	uint32
};

struct SimulationConfig
{
	ESimScalar scalarType = ESimScalar::floatType;
	ESimCount countType = ESimCount::uint32;

	int verbosity;

	std::size_t maxIter;
	std::string gpuSpec;
};

struct Simulation
{
	virtual ~Simulation() = 0;

	virtual void run( SimHostRNG& ) = 0;
	virtual void write_results( input::Parameters const& ) = 0;
};

/** Create `Simulation` instance
 *
 * Creates a `Simulation` instance for the specified `input::Parameters` and
 * `SimulationConfig`.
 */
std::unique_ptr<Simulation> make_simulation(
	SimHostRNG&,
	input::Parameters const&,
	SimulationConfig const& = SimulationConfig{}
);

#endif //  SIMULATION_HPP_37C3055A_87E4_48BB_B6B8_1A4D16482DFC
