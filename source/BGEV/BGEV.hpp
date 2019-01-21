#ifndef BGEV_HPP_
#define BGEV_HPP_

#include <cinttypes>
#include <string>

struct BGEVParameters
{
	uint32_t particleCount;
	int32_t nMax;
	double gamma;
	std::string output;
};

#endif