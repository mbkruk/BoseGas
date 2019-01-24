#include <cstdlib>
#include <immintrin.h>
#include <vector>
#include <iostream>

#include "BGEV.hpp"

__m256d* BGEvolution::derivative(const __m256d *r) 
{
	const int_fast32_t baab_length = 2*nMax+1;
	__m256d baabtab[baab_length];

	//filling baab array
	double a = (*(r))[0];
	double b = (*(r))[1];
	baabtab[nMax] =  _mm256_set_pd(b,a,a,b);

	double c, d;
	for (int_fast32_t i=1;i<=nMax;++i)
	{
		a = (*(r+i))[0];
		b = (*(r+i))[1];
		c = (*(r+i))[2];
		d = (*(r+i))[3];
		baabtab[i+nMax] =  _mm256_set_pd(b,a,a,b);
		baabtab[nMax-i] =  _mm256_set_pd(d,c,c,d);
	}

	//A,B,C,D are used to calculate data for a_k, b_k, a_-k, b_-k
	__m256d A = _mm256_setzero_pd();
	__m256d B = _mm256_setzero_pd();
	__m256d C = _mm256_setzero_pd();
	__m256d D = _mm256_setzero_pd();
	__m256d sgn = _mm256_set_pd(-1.0,1.0,1.0,1.0);

	for (int_fast32_t j=0;j<indices[nMax].size();j+=3) //calculate sum for a_0 and b_0
		{
			A = A+baabtab[indices[nMax][j]+nMax]*_mm256_permute_pd(baabtab[indices[nMax][j+1]+nMax],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax][j+2]+nMax],0b0000)*sgn;
			B = B+_mm256_permute_pd(baabtab[indices[nMax][j]+nMax],0b0101)*_mm256_permute_pd(baabtab[indices[nMax][j+1]+nMax],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax][j+2]+nMax],0b1111)*sgn;
		}

	a = 2.0*gamma*(A[0]+A[1]+A[2]+A[3]);
	b = -2.0*gamma*(B[0]+B[1]+B[2]+B[3]);
	*(pDerivative) = _mm256_set_pd(b,a,b,a);

	for (int_fast32_t i=1;i<=nMax;++i) //calculate  a_i b_i a_-i, b_i 
	{
		A = _mm256_setzero_pd();
		B = _mm256_setzero_pd();
		C = _mm256_setzero_pd();
		D = _mm256_setzero_pd();
		for (int_fast32_t j=0;j<indices[nMax+i].size();j+=3) //calculate sum for a_i and b_i
		{
			A = A+baabtab[indices[nMax+i][j]+nMax]*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]+nMax],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]+nMax],0b0000)*sgn;
			B = B+_mm256_permute_pd(baabtab[indices[nMax+i][j]+nMax],0b0101)*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]+nMax],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]+nMax],0b1111)*sgn;
		}

		for (int_fast32_t j=0;j<indices[nMax-i].size();j+=3) //calculate sum for a_i and b_i
		{
			C = C+baabtab[indices[nMax-i][j]+nMax]*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]+nMax],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]+nMax],0b0000)*sgn;
			D = D+_mm256_permute_pd(baabtab[indices[nMax-i][j]+nMax],0b0101)*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]+nMax],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]+nMax],0b1111)*sgn;
		}

		*(pDerivative+i) = _mm256_set_pd((*(r+i))[2]*i*i*(-1.0)-2.0*gamma*(D[0]+D[1]+D[2]+D[3]),(*(r+i))[3]*i*i+2.0*gamma*(C[0]+C[1]+C[2]+C[3]),
		(*(r+i))[0]*i*i*(-1.0)-2.0*gamma*(B[0]+B[1]+B[2]+B[3]),(*(r+i))[1]*i*i+2.0*gamma*(A[0]+A[1]+A[2]+A[3]));
	}
	return pDerivative;
}

void BGEvolution::rk5()
{
	__m256d k1[stride], k2[stride], k3[stride], k4[stride], k5[stride], k6[stride];
	__m256d s[stride];
	
	for (int_fast32_t i=0;i<=nMax;++i)
		k1[i] = *(derivative(pCurrent)+i)*H;
	for (int_fast32_t i=0;i<=nMax;++i)
		s[i] = *(pCurrent+i)+*(k1+i)*b[1];
 
	for (int_fast32_t i=0;i<=nMax;++i)
		k2[i] = *(derivative(s)+i)*H;
	for (int_fast32_t i=0;i<=nMax;++i)
		s[i] = *(pCurrent+i)+*(k1+i)*b[2]+*(k2+i)*b[3];

	for (int_fast32_t i=0;i<=nMax;++i)
		k3[i] = *(derivative(s)+i)*H;
	for (int_fast32_t i=0;i<=nMax;++i)
		s[i] = *(pCurrent+i)+*(k1+i)*b[4]+*(k2+i)*b[5]+*(k3+i)*b[6];
	
	for (int_fast32_t i=0;i<=nMax;++i)
		k4[i] = *(derivative(s)+i)*H;
	for (int_fast32_t i=0;i<=nMax;++i)
		s[i] = *(pCurrent+i)+*(k1+i)*b[7]+*(k2+i)*b[8]+*(k3+i)*b[9]+*(k4+i)*b[10];

	for (int_fast32_t i=0;i<=nMax;++i)
		k5[i] = *(derivative(s)+i)*H;
	for (int_fast32_t i=0;i<=nMax;++i)
		s[i] = *(pCurrent+i)+*(k1+i)*b[11]+*(k2+i)*b[12]+*(k3+i)*b[13]+*(k4+i)*b[14]+*(k5+i)*b[15];

	for (int_fast32_t i=0;i<=nMax;++i)
		k6[i] = *(derivative(s)+i)*H;

	for (int_fast32_t i=0;i<=nMax;++i)
		*(pCurrent+stride+i) = *(pCurrent+i)+*(k1+i)*c[1]+*(k3+i)*c[3]+*(k4+i)*c[4]+*(k6+i)*c[6];
	
}

void BGEvolution::evolve(uint_fast32_t steps)
{
	while(steps--)
	{
		rk5();
		pCurrent += stride;
	}
	for (int_fast32_t i=0;i<stride;++i)
		initial_conditions[i] = *(pCurrent+i);
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

	H = _mm256_set1_pd(h_);

	particleCount = params.particleCount;
	nMax = params.nMax;
	batchSize = params.batchSize;
	stride = nMax+1;
	gamma = params.gamma;
	
	std::vector<int_fast32_t> A;
	for (int_fast32_t i=-nMax;i<=nMax;++i)
	{
		for (int_fast32_t j=-nMax;j<=nMax;++j)
		for (int_fast32_t k=-nMax;k<=nMax;++k)
		for (int_fast32_t l=-nMax;l<=nMax;++l)
			if (k+l-j==i)
			{
				A.push_back(j);
				A.push_back(k);
				A.push_back(l);
			}
		indices.push_back(A);
		A.clear();
	}
	
	pData = (__m256d*)aligned_alloc(sizeof(__m256d),stride*sizeof(__m256d)*batchSize);
	pCurrent = pData;
	pDerivative = (__m256d*)aligned_alloc(sizeof(__m256d),stride);
	initial_conditions = (__m256d*)aligned_alloc(sizeof(__m256d),stride);
}

void BGEvolution::destroy()
{
	free(pData);
	free(pDerivative);
}

void BGEvolution::stdinICInit()
{
	double a, b;
	for(int_fast32_t i=-nMax;i<=nMax;++i)
	{
		std::cin >> a;
		if (i>0)
			pData[i][0] = a;
		else
		if (i<0)
			pData[abs(i)][2] = a;
		else
		if (i==0)
		{
			pData[0][0] = a;
			pData[0][2] = a;
		}
	}

	for (int_fast32_t i=-nMax;i<=nMax;++i)
	{
		std::cin >> b;
		if (i>0)
			pData[i][1] = b;
		else
		if (i<0)
			pData[abs(i)][3] = b;
		else
		if(i==0)
		{
			pData[0][1] = b;
			pData[0][3] = b;
		}
	}
}