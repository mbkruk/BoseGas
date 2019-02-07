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
	double gamma, h;
	uint_fast32_t batchSize, batchCount;
	std::string output;
};

class BGEvolution
{
private:

	uint_fast32_t particleCount;
	int_fast32_t nMax;
	uint_fast32_t batchSize;
	uint_fast32_t batchCount;
	size_t stride;
	double gamma;
	double h;
	std::string output;
	std::string interactionType;
	std::vector<double> interactionCoefficients;

	// non-zero interaction sum factors
	std::vector<std::vector<int_fast32_t> > indices;

	//initial conditions
	__m256d *initial_conditions; 
	__m256d *baabtab;

	// RK5 parameters
	__m256d b[16];
	__m256d c[7];
	__m256d H;

	// {{re_0,im_0,re_0,im_0}, {re_1,im_1,re_-1,im_-1}, ...}
	__m256d *pData;
	__m256d *pCurrent;
	__m256d *pDerivative;

	__m256d *k1, *k2, *k3, *k4, *k5, *k6;
	__m256d *yk1, *yk2, *yk3, *yk4, *yk5;

	void derivative(const __m256d *r);

public:
	void evolve(uint_fast32_t steps);

	void create(const BGEVParameters &params);
	void destroy();
	void icInit();
	
	void stdinICInit();
	void printParameters();

	double averageNZero();
	double fluctuationsNZero(const double n0);

	double nAll();
	double nZero();
	double momentum();
	double kineticEnergy();
	double potentialEnergy();
};

#endif
