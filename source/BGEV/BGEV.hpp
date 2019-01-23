#ifndef BGEV_HPP_
#define BGEV_HPP_

#include <cinttypes>
#include <string>
#include <vector>
#include <immintrin.h>

struct BGEVParameters
{
	uint_fast32_t particleCount;
	int_fast32_t nMax;
	double gamma;
	uint_fast32_t batchSize;
	std::string output;
};

class BGEvolution
{
private:

	uint_fast32_t particleCount;
	int_fast32_t nMax;
	uint_fast32_t batchSize;
	size_t stride;
	double gamma;

	// non-zero interaction sum factors
	std::vector<std::vector<int_fast32_t> > indices;

	// RK5 parameters
	__m256d b[16];
	__m256d c[7];
	__m256d h;

	// {{re_0,im_0,re_0,im_0}, {re_1,im_1,re_-1,im_-1}, ...}
	__m256d *pData;
	__m256d *pCurrent;
	__m256d *pDerivative;

public:

	void evolve(uint_fast32_t steps);

	void create(double h_, const BGEVParameters &params);
	void destroy();
};

#endif