#ifndef BGMC_HPP_
#define BGMC_HPP_

#include <cinttypes>
#include <string>

struct BGMCParameters
{
	uint32_t particleCount;
	int32_t nMax;
	double beta;
	std::string output;
};

#endif
