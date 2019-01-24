#include <cstdlib>
#include <immintrin.h>
#include <cmath>
#include <iostream>

#include "BGEV.hpp"

double BGEvolution::averageNZero(const int_fast32_t steps)
{	
	double aveN0;
	int_fast32_t int_steps;

	if (steps%2==0) //we want even number of steps
		int_steps = steps;
	else
		int_steps = steps-1;

	int_fast32_t avx_steps = int_steps/2-1;

	double integral = (*(pData))[0]*(*(pData))[0]+(*(pData))[1]*(*(pData))[1]; //initialization with N0 value at t=0
	__m256d avx_sum = _mm256_set1_pd(0.0);
	__m256d avx_two = _mm256_set1_pd(2.0);
	for (int_fast32_t i=1;i<=avx_steps;++i)
		avx_sum += (*(pData+2*i*stride))*(*(pData+2*i*stride))+avx_two*(*(pData+(2*i-1)*stride))*(*(pData+(2*i-1)*stride));
	avx_sum += avx_two*(*(pData+(int_steps-1)*stride))*(*(pData+(int_steps-1)*stride));
	integral += avx_sum[0]+avx_sum[1]+avx_sum[2]+avx_sum[3]+(*(pData+int_steps*stride))[0]*(*(pData+int_steps*stride))[0]+(*(pData+int_steps*stride))[1]*(*(pData+int_steps*stride))[1];
	aveN0 = integral/(int_steps*3.0);

	for(int i=0;i<nMax+1;++i)
	{
		std::cerr << (*(pData+i))[0] << '\n';
	}
	return aveN0;
}

double BGEvolution::fluctuationsNZero(const int_fast32_t steps, const double n0)
{
	double aveN0sq;
	int_fast32_t int_steps;

	if (steps%2==0) //we want even number of steps
		int_steps = steps;
	else
		int_steps = steps-1;

	int_fast32_t avx_steps = int_steps/2-1;

	double integral = pow((*(pData))[0],4)+pow((*(pData))[1],4)+2.0*pow((*(pData))[0],2)*pow((*(pData))[1],2); //initialization with N0^2 value at t=0.
	__m256d avx_sum = _mm256_set1_pd(0.0);
	__m256d avx_half = _mm256_set1_pd(0.5);
	__m256d avx_two = _mm256_set1_pd(2.0);

	__m256d X1,X2,Y1,Y2;
	for (int_fast32_t i=1;i<=avx_steps;++i)
	{
		X1 = (*(pData+2*i*stride))*(*(pData+2*i*stride)); 
		X2 = _mm256_permute_pd(X1, 0b0101);
		Y1 = (*(pData+(2*i-1)*stride))*(*(pData+(2*i-1)*stride));
		Y2 = _mm256_permute_pd(Y1, 0b0101);
	
		avx_sum += avx_half*(X1*X1+avx_two*X1*X2+X2*X2)+Y1*Y1+avx_two*Y1*Y2+Y2*Y2;
	}	

	Y1 = (*(pData+(int_steps-1)*stride))*(*(pData+(int_steps-1)*stride));
	Y2 = _mm256_permute_pd(Y1, 0b0101);

	avx_sum += Y1*Y1+avx_two*Y1*Y2+Y2*Y2;
	integral += avx_sum[0]+avx_sum[1]+avx_sum[2]+avx_sum[3]+pow((*(pData+int_steps*stride))[0],4)+pow((*(pData+int_steps*stride))[1],4)+
	2.0*pow((*(pData+int_steps*stride))[0],2)*pow((*(pData+int_steps*stride))[1],2);
	
	integral *= h/3.0;

	aveN0sq = integral/(h*int_steps);
	return sqrt(aveN0sq-n0*n0);
}

