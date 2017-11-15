#ifndef OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D
#define OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D

#include <string>
#include <vector>

#include <cstddef>

#include "common.hpp"

namespace output
{
	struct Output
	{
		EModel model;

		std::size_t componentCount;
		std::size_t abcSampleCount;

		std::vector<double> m;
		std::vector<double> s;
		std::vector<double> c;
		std::vector<double> az;
		std::vector<double> dist;
		std::vector<double> w;

		double epsilon;

		// optional meta data
		std::string producerName;
		std::string date;
		std::string machine;
	};

	Output load( char const* );
}

#endif // OUTPUT_HPP_951E9419_17C4_432F_AAA5_FB7D1A153C1D
