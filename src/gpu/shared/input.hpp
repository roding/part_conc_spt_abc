#ifndef INPUT_HPP_C2353B1D_0C7B_4B27_832E_F970C522DBCB
#define INPUT_HPP_C2353B1D_0C7B_4B27_832E_F970C522DBCB

#include <string>
#include <vector>

#include <cstddef>

#include "common.hpp"

namespace input
{
	struct Bounds
	{
		double lower;
		double upper;
	};
	
	struct Parameters
	{
		EModel model;
		EWeightingScheme weightingScheme;

		std::size_t componentCount;
		std::size_t kmin;
		std::size_t deBinCount;
		std::size_t abcSampleCount;
		std::size_t avgTrialCount;

		double Lx, Ly, Lz;
		double gamma;
		double deltaGamma;
		bool adaptiveGamma;

		Bounds de, m, s, c, az;
		
		std::string outputFilePath;

		
		double ax, ay;
		double deltaT;
		std::vector<std::size_t> frameCounts;
		std::vector<std::size_t> Ks;
		std::vector<double> DEs;
	};

	/** Load parameters from XML input
	 *
	 * Loads the simulation parameters from the specified XML; this should
	 * correspond to the `input.xml` compared to the Julia code. The optional
	 * second argument may be used to override the value of `data_file_path`,
	 * causing the function to load a different `data.xml`.
	 */
	Parameters load( char const*, char const* = nullptr );
}

#endif // INPUT_HPP_C2353B1D_0C7B_4B27_832E_F970C522DBCB
