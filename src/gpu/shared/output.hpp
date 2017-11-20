#ifndef OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D
#define OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D

#include <string>
#include <vector>

#include <unordered_set>
#include <unordered_map>

#include <cstddef>

#include "common.hpp"

namespace output
{
	struct Output
	{
		EModel model;
		EWeightingScheme weightingScheme;

		std::size_t abcSampleCount;
		std::size_t componentCount, zComponentCount;

		std::vector<double> m;
		std::vector<double> c;
		std::vector<double> az;
		std::vector<double> dist;
		std::vector<double> w;

		bool converged;
		double epsilon;

		std::unordered_map<std::string,std::string> meta;
	};

	void write( char const*, Output const& );

	Output load( char const* ); // XXX-deprecated, replace with below

	/* TODO:
	Output load( 
		char const*, 
		std::unordered_set<std::string> const& = { "app", "date", "machine" }
	);
	*/
}

#endif // OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D
