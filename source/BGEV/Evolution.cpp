#include <cstdlib>

#include "BGEV.hpp"

void BGEvolution::evolve(uint_fast32_t steps)
{
	while(steps--)
	{
		pCurrent += stride;
	}
}

void BGEvolution::create(double h_, const BGEVParameters &params)
{
	b[1] = _mm256_set1_pd(0.2);
	b[2] = _mm256_set1_pd(3.0/40.0);
	b[3] = _mm256_set1_pd(9.0/40.0);
	b[4] = _mm256_set1_pd(0.3);
	b[5] = _mm256_set1_pd(-0.9);
	b[6] = _mm256_set1_pd(6.0/5.0);
	b[7] = _mm256_set1_pd(-11.0/54.0);
	b[8] = _mm256_set1_pd(2.5);
	b[9] = _mm256_set1_pd(-70.0/27.0);
	b[10] = _mm256_set1_pd(35.0/27.0);
	b[11] = _mm256_set1_pd(1631.0/55296.0);
	b[12] = _mm256_set1_pd(175.0/512.0);
	b[13] = _mm256_set1_pd(575.0/13824.0);
	b[14] = _mm256_set1_pd(44275.0/110592.0);
	b[15] = _mm256_set1_pd(253.0/4096.0);

	c[1] = _mm256_set1_pd(37.0/378.0);
	c[2] = _mm256_set1_pd(0.0);
	c[3] = _mm256_set1_pd(250.0/621.0);
	c[4] = _mm256_set1_pd(125.0/594.0);
	c[5] = _mm256_set1_pd(0.0);
	c[6] = _mm256_set1_pd(512.0/1771.0);

	h = _mm256_set1_pd(h_);

	particleCount = params.particleCount;
	nMax = params.nMax;
	batchSize = params.batchSize;
	stride = nMax+1;
	gamma = params.gamma;

	pData = (__m256d*)aligned_alloc(sizeof(__m256d),stride*sizeof(__m256d)*batchSize);
	pCurrent = pData;
	pDerivative = (__m256d*)aligned_alloc(sizeof(__m256d),stride);
}

void BGEvolution::destroy()
{
	free(pData);
	free(pDerivative);
}
