#include <immintrin.h>
#include <cmath>

#include "BGEV.hpp"

double BGEvolution::nAll()
{
	double n = (*(pCurrent))[0]*(*(pCurrent))[0]+(*(pCurrent))[1]*(*(pCurrent))[1];
	for (int_fast32_t i=1;i<=nMax;++i)
		n += (*(pCurrent+i))[0]*(*(pCurrent+i))[0]+(*(pCurrent+i))[1]*(*(pCurrent+i))[1]+(*(pCurrent+i))[2]*(*(pCurrent+i))[2]+(*(pCurrent+i))[3]*(*(pCurrent+i))[3]; 

	return n;
}

double BGEvolution::nZero()
{
	return (*(pCurrent))[0]*(*(pCurrent))[0]+(*(pCurrent))[1]*(*(pCurrent))[1];
}

double BGEvolution::momentum()
{
	double p = 0.0;
	for (int_fast32_t i=1;i<=nMax;++i)
		p += i*((*(pCurrent+i))[0]*(*(pCurrent+i))[0]+(*(pCurrent+i))[1]*(*(pCurrent+i))[1])-i*((*(pCurrent+i))[2]*(*(pCurrent+i))[2]+(*(pCurrent+i))[3]*(*(pCurrent+i))[3]); 

	return p;
}

double BGEvolution::kineticEnergy()
{
	double e = 0.0;
	for (int_fast32_t i=1;i<=nMax;++i)
		e += i*i*((*(pCurrent+i))[0]*(*(pCurrent+i))[0]+(*(pCurrent+i))[1]*(*(pCurrent+i))[1]+(*(pCurrent+i))[2]*(*(pCurrent+i))[2]+(*(pCurrent+i))[3]*(*(pCurrent+i))[3]);

	return e;
}

int sgn(int x)
{
	if (x>=0)
		return 1;
	return -1;
}

double BGEvolution::potentialEnergy()
{
	double u = 0.0;
	if (interactionType=="contact")
	{
	for (int_fast32_t i=-nMax;i<=nMax;++i)
	for (int_fast32_t j=-nMax;j<=nMax;++j)
	for (int_fast32_t k=-nMax;k<=nMax;++k)
	for (int_fast32_t l=-nMax;l<=nMax;++l)	
		if(i+j==k+l)
			u += (*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]-
			(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
			(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
			(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]+
			(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
			(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]-
			(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]+
			(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)];
		
	return gamma*u;
	}
	else
	{
		for (int_fast32_t i=-nMax;i<=nMax;++i)
		for (int_fast32_t j=-nMax;j<=nMax;++j)
		for (int_fast32_t k=-nMax;k<=nMax;++k)
		for (int_fast32_t l=-nMax;l<=nMax;++l)	
			if(i+j==k+l)
				u += reducedCoefficients[std::max(abs(k-i),abs(l-j))]*((*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]-
				(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
				(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
				(*(pCurrent+abs(i)))[1-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]+
				(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]+
				(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[1-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]-
				(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[1-sgn(k)]*(*(pCurrent+abs(l)))[1-sgn(l)]+
				(*(pCurrent+abs(i)))[2-sgn(i)]*(*(pCurrent+abs(j)))[2-sgn(j)]*(*(pCurrent+abs(k)))[2-sgn(k)]*(*(pCurrent+abs(l)))[2-sgn(l)]);
			
		return gamma*u;
	}


}