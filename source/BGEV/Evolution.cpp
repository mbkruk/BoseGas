#include <cstdlib>
#include <immintrin.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "BGEV.hpp"

const int_fast32_t dist = 2*sizeof(double)+9;

#ifndef _ISOC11_SOURCE
#ifdef _WIN32
#define aligned_alloc(alignment,size) _aligned_malloc((size),(alignment))
#elif _POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600
static void* aligned_alloc(size_t alignment, size_t size)
{
	void *ptr;
	posix_memalign(&ptr,alignment,size);
	return ptr;
}
#endif
#endif

void BGEvolution::derivative(const __m256d *r)
{
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
			A = A+baabtab[indices[nMax][j]]*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b0000)*sgn;
			B = B+_mm256_permute_pd(baabtab[indices[nMax][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b1111)*sgn;
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
			A = A+baabtab[indices[nMax+i][j]]*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]],0b0000)*sgn;
			B = B+_mm256_permute_pd(baabtab[indices[nMax+i][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]],0b1111)*sgn;
		}

		for (int_fast32_t j=0;j<indices[nMax-i].size();j+=3) //calculate sum for a_-i and b_-i
		{
			C = C+baabtab[indices[nMax-i][j]]*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]],0b0000)*sgn;
			D = D+_mm256_permute_pd(baabtab[indices[nMax-i][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]],0b1111)*sgn;
		}

		*(pDerivative+i) = _mm256_set_pd((*(r+i))[2]*i*i*(-1.0)-2.0*gamma*(D[0]+D[1]+D[2]+D[3]),(*(r+i))[3]*i*i+2.0*gamma*(C[0]+C[1]+C[2]+C[3]),
		(*(r+i))[0]*i*i*(-1.0)-2.0*gamma*(B[0]+B[1]+B[2]+B[3]),(*(r+i))[1]*i*i+2.0*gamma*(A[0]+A[1]+A[2]+A[3]));
	}
}

void BGEvolution::derivativeLong(const __m256d *r)
{
	//derivative with presence of long-interaction non-trivial coefficients
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

	__m256d A = _mm256_setzero_pd();
	__m256d B = _mm256_setzero_pd();
	__m256d C = _mm256_setzero_pd();
	__m256d D = _mm256_setzero_pd();
	__m256d sgn = _mm256_set_pd(-1.0,1.0,1.0,1.0);
	for (int_fast32_t j=0;j<indices[nMax].size();j+=3)
	{
		A = A+interactionCoefficients[indicesCount[nMax]+j/3]*baabtab[indices[nMax][j]]*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b0110)*
		_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b0000)*sgn;
		B = B+interactionCoefficients[indicesCount[nMax]+j/3]*_mm256_permute_pd(baabtab[indices[nMax][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b1001)*
		_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b1111)*sgn;
	}
	a = 2.0*gamma*(A[0]+A[1]+A[2]+A[3]);
	b = -2.0*gamma*(B[0]+B[1]+B[2]+B[3]);
	*(pDerivative) = _mm256_set_pd(b,a,b,a);

	for (int_fast32_t i=1;i<=nMax;++i)
	{
		A = _mm256_setzero_pd();
		B = _mm256_setzero_pd();
		C = _mm256_setzero_pd();
		D = _mm256_setzero_pd();
		for (int_fast32_t j=0;j<indices[nMax+i].size();j+=3)
		{
			A = A+interactionCoefficients[indicesCount[nMax+i]+j/3]*baabtab[indices[nMax+i][j]]*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]],0b0000)*sgn;
			B = B+interactionCoefficients[indicesCount[nMax+i]+j/3]*_mm256_permute_pd(baabtab[indices[nMax+i][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax+i][j+1]],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax+i][j+2]],0b1111)*sgn;
		}

		for (int_fast32_t j=0;j<indices[nMax-i].size();j+=3)
		{
			C = C+interactionCoefficients[indicesCount[nMax-i]+j/3]*baabtab[indices[nMax-i][j]]*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]],0b0110)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]],0b0000)*sgn;
			D = D+interactionCoefficients[indicesCount[nMax-i]+j/3]*_mm256_permute_pd(baabtab[indices[nMax-i][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax-i][j+1]],0b1001)*
			_mm256_permute_pd(baabtab[indices[nMax-i][j+2]],0b1111)*sgn;
		}

		*(pDerivative+i) = _mm256_set_pd((*(r+i))[2]*i*i*(-1.0)-2.0*gamma*(D[0]+D[1]+D[2]+D[3]),(*(r+i))[3]*i*i+2.0*gamma*(C[0]+C[1]+C[2]+C[3]),
		(*(r+i))[0]*i*i*(-1.0)-2.0*gamma*(B[0]+B[1]+B[2]+B[3]),(*(r+i))[1]*i*i+2.0*gamma*(A[0]+A[1]+A[2]+A[3]));
	}
}

void BGEvolution::evolve()
{
	uint_fast32_t steps = batchSize;
	if (interactionType=="custom")
	{

		while (steps--)
		{
			derivativeLong(pCurrent);
			for (int_fast32_t i=0;i<=nMax;++i)
				k1[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk1[i] = *(pCurrent+i)+*(k1+i)*b[1];

			derivativeLong(yk1);
			for (int_fast32_t i=0;i<=nMax;++i)
				k2[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk2[i] = *(pCurrent+i)+*(k1+i)*b[2]+*(k2+i)*b[3];

			derivativeLong(yk2);
			for (int_fast32_t i=0;i<=nMax;++i)
				k3[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk3[i] = *(pCurrent+i)+*(k1+i)*b[4]+*(k2+i)*b[5]+*(k3+i)*b[6];

			derivativeLong(yk3);
			for (int_fast32_t i=0;i<=nMax;++i)
				k4[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk4[i] = *(pCurrent+i)+*(k1+i)*b[7]+*(k2+i)*b[8]+*(k3+i)*b[9]+*(k4+i)*b[10];

			derivativeLong(yk4);
			for (int_fast32_t i=0;i<=nMax;++i)
				k5[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk5[i] = *(pCurrent+i)+*(k1+i)*b[11]+*(k2+i)*b[12]+*(k3+i)*b[13]+*(k4+i)*b[14]+*(k5+i)*b[15];

			derivativeLong(yk5);
			for (int_fast32_t i=0;i<=nMax;++i)
				k6[i] = *(pDerivative+i)*H;

			for (int_fast32_t i=0;i<=nMax;++i)
				*(pCurrent+stride+i) = *(pCurrent+i)+*(k1+i)*c[1]+*(k3+i)*c[3]+*(k4+i)*c[4]+*(k6+i)*c[6];

			pCurrent += stride;
		}
	}
	else
	if (interactionType=="contact")
	{
		while (steps--)
		{
			derivative(pCurrent);
			for (int_fast32_t i=0;i<=nMax;++i)
				k1[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk1[i] = *(pCurrent+i)+*(k1+i)*b[1];

			derivative(yk1);
			for (int_fast32_t i=0;i<=nMax;++i)
				k2[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk2[i] = *(pCurrent+i)+*(k1+i)*b[2]+*(k2+i)*b[3];

			derivative(yk2);
			for (int_fast32_t i=0;i<=nMax;++i)
				k3[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk3[i] = *(pCurrent+i)+*(k1+i)*b[4]+*(k2+i)*b[5]+*(k3+i)*b[6];

			derivative(yk3);
			for (int_fast32_t i=0;i<=nMax;++i)
				k4[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk4[i] = *(pCurrent+i)+*(k1+i)*b[7]+*(k2+i)*b[8]+*(k3+i)*b[9]+*(k4+i)*b[10];

			derivative(yk4);
			for (int_fast32_t i=0;i<=nMax;++i)
				k5[i] = *(pDerivative+i)*H;
			for (int_fast32_t i=0;i<=nMax;++i)
				yk5[i] = *(pCurrent+i)+*(k1+i)*b[11]+*(k2+i)*b[12]+*(k3+i)*b[13]+*(k4+i)*b[14]+*(k5+i)*b[15];

			derivative(yk5);
			for (int_fast32_t i=0;i<=nMax;++i)
				k6[i] = *(pDerivative+i)*H;

			for (int_fast32_t i=0;i<=nMax;++i)
				*(pCurrent+stride+i) = *(pCurrent+i)+*(k1+i)*c[1]+*(k3+i)*c[3]+*(k4+i)*c[4]+*(k6+i)*c[6];

			pCurrent += stride;
		}
	}

	for (int_fast32_t i=0;i<stride;++i)
		*(initial_conditions+i) = *(pCurrent+i);
}

void BGEvolution::create(const BGEVParameters &params)
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

	H = _mm256_set1_pd(params.h);
	h = params.h;

	particleCount = params.particleCount;
	nMax = params.nMax;
	batchSize = params.batchSize;
	batchCount = params.batchCount;
	alphaCount = params.alphaCount;
	stride = nMax+1;
	gamma = params.gamma;
	output = params.output;
	interactionType = params.interactionType;

	for (int_fast32_t i=0;i<params.reducedCoefficients.size();++i)
		reducedCoefficients.push_back(params.reducedCoefficients[i]);

	std::vector<int_fast32_t> A;
	if (interactionType=="custom")
	{
		indicesCount = (int_fast32_t*)aligned_alloc(sizeof(int_fast32_t),(2*nMax+2)*sizeof(int_fast32_t));
		indicesCount[0] = 0;

		for (int_fast32_t i=-nMax;i<=nMax;++i)
		{
			indicesCount[i+1+nMax] = indicesCount[i+nMax];
			for (int_fast32_t j=-nMax;j<=nMax;++j)
			for (int_fast32_t k=-nMax;k<=nMax;++k)
			for (int_fast32_t l=-nMax;l<=nMax;++l)
				if (k+l-j==i)
				{
					A.push_back(j+nMax);
					A.push_back(k+nMax);
					A.push_back(l+nMax);

					++indicesCount[i+1+nMax];
				}

			indices.push_back(A);
			A.clear();
		}

		interactionCoefficients = (__m256d*)aligned_alloc(sizeof(__m256d),(indicesCount[2*nMax+1])*sizeof(__m256d));
		__m256d *p;
		p = interactionCoefficients;
		for (int_fast32_t i=-nMax;i<=nMax;++i)
		for (int_fast32_t j=-nMax;j<=nMax;++j)
		for (int_fast32_t k=-nMax;k<=nMax;++k)
		for (int_fast32_t l=-nMax;l<=nMax;++l)
			if (k+l-j==i)
				{
					*p = _mm256_set1_pd(reducedCoefficients[abs(k-j)]);
					p++;
				}
	}
	else
	if (interactionType=="contact")
	{
		for (int_fast32_t i=-nMax;i<=nMax;++i)
		{
			for (int_fast32_t j=-nMax;j<=nMax;++j)
			for (int_fast32_t k=-nMax;k<=nMax;++k)
			for (int_fast32_t l=-nMax;l<=nMax;++l)
				if (k+l-j==i)
				{
					A.push_back(j+nMax);
					A.push_back(k+nMax);
					A.push_back(l+nMax);
				}

			indices.push_back(A);
			A.clear();
		}
	}

	pData = (__m256d*)aligned_alloc(sizeof(__m256d),stride*sizeof(__m256d)*(batchSize+1));
	pDerivative = (__m256d*)aligned_alloc(sizeof(__m256d),stride*sizeof(__m256d));
	k1 = (__m256d*)aligned_alloc(sizeof(__m256d),11*stride*sizeof(__m256d));
	initial_conditions = (__m256d*)aligned_alloc(sizeof(__m256d),stride*sizeof(__m256d));
	baabtab = (__m256d*)aligned_alloc(sizeof(__m256d),(2*nMax+1)*sizeof(__m256d));

	k2 = k1+stride;
	k3 = k2+stride;
	k4 = k3+stride;
	k5 = k4+stride;
	k6 = k5+stride;
	yk1 = k6+stride;
	yk2 = yk1+stride;
	yk3 = yk2+stride;
	yk4 = yk3+stride;
	yk5 = yk4+stride;

	threadCount = params.threadCount;

	if (threadCount>0)
	{
		pThreads = new std::thread[threadCount];

		barrier.setThreadCount(threadCount+1);

		for (uint32_t i=0;i<threadCount;++i)
			pThreads[i] = std::thread(&BGEvolution::thread,this,i+1);
	}
	else
		pThreads = nullptr;
}

void BGEvolution::destroy()
{
	if (pThreads)
	{
		for (uint32_t i=0;i<threadCount;++i)
			pThreads[i].join();
		delete [] pThreads;
		pThreads = nullptr;
		threadCount = 0;
	}
	free(pData);
	free(pDerivative);
	free(k1);
}

void BGEvolution::stdinInit()
{
	double a, b;
	for(int_fast32_t i=-nMax;i<=nMax;++i)
	{
		std::cin >> a;
		if (i>0)
			(*(initial_conditions+i))[0] = a;
		else
		if (i<0)
			(*(initial_conditions+abs(i)))[2] = a;
		else
		if (i==0)
		{
			(*(initial_conditions))[0] = a;
			(*(initial_conditions))[2] = a;
		}
	}

	for (int_fast32_t i=-nMax;i<=nMax;++i)
	{
		std::cin >> b;
		if (i>0)
			(*(initial_conditions+i))[1] = b;
		else
		if (i<0)
			(*(initial_conditions+abs(i)))[3] = b;
		else
		if (i==0)
		{
			(*(initial_conditions))[1] = b;
			(*(initial_conditions))[3] = b;
		}
	}

	std::fstream output_file;
	output_file.open(output,std::ios::out);

	output_file << "Number of particles: " << particleCount << '\n';
	output_file << "Nmax: " << nMax << '\n';
	output_file << "Gamma: " << gamma << '\n';
	output_file << "Interaction type: " << interactionType << '\n';
	output_file << "Step size: " << h << '\n';
	output_file << "Batch size: " << batchSize << '\n';
	output_file << "Batch count: " << batchCount << '\n';
	output_file << std::setw(dist) << "Batch number" << std::setw(dist) << std::setprecision(10) << "Average n0"
		<< std::setw(dist) << "Fluctuations n0" << std::setw(dist) << "Energy"
		<<  std::setw(dist) << "Particle count" << std::setw(dist) << "Momentum";
	if(interactionType=="contact")
		output_file << std::setw(dist) << "Constant X" << '\n';
	else
		output_file << '\n';
	output_file.close();
}

void BGEvolution::icInit()
{
	for (int_fast32_t i=0;i<stride;++i)
		*(pData+i) = *(initial_conditions+i);

	pCurrent = pData;
}

void BGEvolution::printParameters()
{
	std::cerr << "Number of particles: " << particleCount << '\n';
	std::cerr << "Nmax: " << nMax << '\n';
	std::cerr << "Gamma: " << gamma << '\n';
	std::cerr << "Interaction type: " << interactionType << '\n';
	std::cerr << "Step size: " << h << '\n';
	std::cerr << "Batch size: " << batchSize << '\n';
	std::cerr << "Batch count: " << batchCount << '\n' << '\n';
	std::cerr << std::setw(dist) << "Batch number" << std::setw(dist) << std::setprecision(10) << "Average n0"
		<< std::setw(dist) << "Fluctuations n0" << std::setw(dist) << "Energy"
		<<  std::setw(dist) << "Particle count" << std::setw(dist) << "Momentum";
	if(interactionType=="contact")
		std::cerr << std::setw(dist) << "Constant X" << '\n';
	else
		std::cerr << '\n';
}

void BGEvolution::saveToFile(const int_fast32_t i)
{
	std::fstream output_file;
	output_file.open(output,std::ios::app|std::ios::out);

	output_file << std::setw(dist) << i+1 << std::setw(dist) << std::setprecision(16) << avgs.back() << std::setw(dist) << flucs.back()
		<< std::setw(dist) << kineticEnergy()+potentialEnergy() <<  std::setw(dist) << nAll() << std::setw(dist) << momentum();
	if(interactionType=="contact")
		output_file << std::setw(dist) << constantx() << '\n';
	else
		output_file << '\n';
	output_file.close();
}

void BGEvolution::lastBatch()
{
	//We have to save last point of evolution (to be able to continue evolving if we want to) and calculate final averages.
	double finalavg, finalfluc;
	double suma = 0.0;
	double sumf = 0.0;
	for (int_fast32_t i=0;i<batchCount;++i)
	{
		suma += avgs[i];
		sumf += flucs[i]*flucs[i]+avgs[i]*avgs[i];
	}

	finalavg = suma/batchCount;
	finalfluc = sqrt(sumf/batchCount-finalavg*finalavg);

	std::cerr << '\n' << "Final average n0: "  << std::setprecision(16) << finalavg << '\n';
	std::cerr << "Final fluctuations n0: " << std::setprecision(16) << finalfluc << '\n';

	std::fstream output_file;
	output_file.open(output,std::ios::app|std::ios::out);

	output_file << '\n' << "Final average n0: "  << std::setprecision(16) << finalavg << '\n';
	output_file << "Final fluctuations n0: " << std::setprecision(16) << finalfluc << '\n';

	output_file << '\n' << "Averages:" << '\n';
	for (int_fast32_t i=0;i<batchCount;++i)
		output_file << avgs[i] << '\n';

	output_file << '\n' << "Fluctuations:" << '\n';
	for (int_fast32_t i=0;i<batchCount;++i)
		output_file << flucs[i] << '\n';

	output_file << '\n' << "Last point and parameters for further simulation:" << '\n';
	output_file << particleCount << '\n';
	output_file << nMax << '\n';
	output_file << gamma << '\n';
	output_file << interactionType << '\n';
	output_file << alphaCount << '\n' << '\n';


	if (interactionType=="custom")
		for (int_fast32_t i=0;i<=2*nMax;++i)
			output_file << reducedCoefficients[i] << '\n';

	output_file << '\n';
	for (int_fast32_t i=-nMax;i<=nMax;++i)
	{
		if (i<0)
			output_file << (*(pCurrent+abs(i)))[2] << '\n';
		else
			output_file << (*(pCurrent+i))[0] << '\n';
	}

	for (int_fast32_t i=-nMax;i<=nMax;++i)
	{
		if (i<0)
			output_file << (*(pCurrent+abs(i)))[3] << '\n';
		else
			output_file << (*(pCurrent+i))[1] << '\n';
	}

	output_file << '\n' << "Final quantities:" << '\n';
	output_file << std::setprecision(16) << finalavg << '\n';
	output_file << std::setprecision(16) << finalfluc << '\n';
	output_file.close();
}
