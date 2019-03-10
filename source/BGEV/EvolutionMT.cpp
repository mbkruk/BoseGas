#include "BGEV.hpp"

void BGEvolution::derivative0(const __m256d *r)
{
	double a = (*(r))[0];
	double b = (*(r))[1];
	baabtab[nMax] =  _mm256_set_pd(b,a,a,b);

	__m256d A = _mm256_setzero_pd();
	__m256d B = _mm256_setzero_pd();
	__m256d C = _mm256_setzero_pd();
	__m256d D = _mm256_setzero_pd();
	__m256d sgn = _mm256_set_pd(-1.0,1.0,1.0,1.0);

	barrier.wait();

	for (int_fast32_t j=0;j<indices[nMax].size();j+=3) //calculate sum for a_0 and b_0
	{
		A = A+baabtab[indices[nMax][j]]*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b0110)*
		_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b0000)*sgn;
		B = B+_mm256_permute_pd(baabtab[indices[nMax][j]],0b0101)*_mm256_permute_pd(baabtab[indices[nMax][j+1]],0b1001)*
		_mm256_permute_pd(baabtab[indices[nMax][j+2]],0b1111)*sgn;
	}

	a = 2.0*gamma*(A[0]+A[1]+A[2]+A[3]);
	b = -2.0*gamma*(B[0]+B[1]+B[2]+B[3]);
	pDerivative[0] = _mm256_set_pd(b,a,b,a);
}

void BGEvolution::derivative1(const __m256d *r, const uint_fast32_t i)
{
	double a, b, c, d;
	a = r[i][0];
	b = r[i][1];
	c = r[i][2];
	d = r[i][3];
	baabtab[i+nMax] =  _mm256_set_pd(b,a,a,b);
	baabtab[nMax-i] =  _mm256_set_pd(d,c,c,d);

	barrier.wait();

	__m256d A = _mm256_setzero_pd();
	__m256d B = _mm256_setzero_pd();
	__m256d C = _mm256_setzero_pd();
	__m256d D = _mm256_setzero_pd();
	__m256d sgn = _mm256_set_pd(-1.0,1.0,1.0,1.0);

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

	pDerivative[i] = _mm256_set_pd(r[i][2]*i*i*(-1.0)-2.0*gamma*(D[0]+D[1]+D[2]+D[3]),r[i][3]*i*i+2.0*gamma*(C[0]+C[1]+C[2]+C[3]),
	r[i][0]*i*i*(-1.0)-2.0*gamma*(B[0]+B[1]+B[2]+B[3]),r[i][1]*i*i+2.0*gamma*(A[0]+A[1]+A[2]+A[3]));
}

void BGEvolution::thread(const uint32_t i)
{
	uint_fast32_t batches = batchCount;
	while (--batches)
	{
		uint_fast32_t steps = batchSize;
		while (--steps)
		{
			derivative1(pCurrent,i);
			k1[i] = pDerivative[i]*H;
			yk1[i] = pCurrent[i]+k1[i]*b[1];

			derivative1(yk1,i);
			k2[i] = *(pDerivative+i)*H;
			yk2[i] = *(pCurrent+i)+*(k1+i)*b[2]+*(k2+i)*b[3];

			derivative1(yk2,i);
			k3[i] = *(pDerivative+i)*H;
			yk3[i] = *(pCurrent+i)+*(k1+i)*b[4]+*(k2+i)*b[5]+*(k3+i)*b[6];

			derivative1(yk3,i);
			k4[i] = *(pDerivative+i)*H;
			yk4[i] = *(pCurrent+i)+*(k1+i)*b[7]+*(k2+i)*b[8]+*(k3+i)*b[9]+*(k4+i)*b[10];

			derivative1(yk4,i);
			k5[i] = *(pDerivative+i)*H;
			yk5[i] = *(pCurrent+i)+*(k1+i)*b[11]+*(k2+i)*b[12]+*(k3+i)*b[13]+*(k4+i)*b[14]+*(k5+i)*b[15];

			derivative1(yk5,i);
			k6[i] = *(pDerivative+i)*H;

			pCurrent[stride+i] = *(pCurrent+i)+*(k1+i)*c[1]+*(k3+i)*c[3]+*(k4+i)*c[4]+*(k6+i)*c[6];
		}
		barrier.wait();
	}
}

void BGEvolution::evolve2()
{
	uint_fast32_t steps = batchSize;
	constexpr uint_fast32_t i = 0;
	while (--steps)
	{
		derivative0(pCurrent);
		k1[0] = pDerivative[0]*H;
		yk1[0] = pCurrent[0]+k1[0]*b[1];

		derivative0(yk1);
		k2[0] = pDerivative[0]*H;
		yk2[0] = pCurrent[0]+k1[0]*b[2]+*(k2+i)*b[3];

		derivative0(yk2);
		k3[0] = pDerivative[0]*H;
		yk3[0] = pCurrent[0]+k1[0]*b[4]+*(k2+i)*b[5]+*(k3+i)*b[6];

		derivative0(yk3);
		k4[0] = pDerivative[0]*H;
		yk4[0] = pCurrent[0]+k1[0]*b[7]+*(k2+i)*b[8]+*(k3+i)*b[9]+*(k4+i)*b[10];

		derivative0(yk4);
		k5[0] = pDerivative[0]*H;
		yk5[0] = pCurrent[0]+k1[0]*b[11]+*(k2+i)*b[12]+*(k3+i)*b[13]+*(k4+i)*b[14]+*(k5+i)*b[15];

		derivative0(yk5);
		k6[0] = pDerivative[0]*H;

		pCurrent[stride+0] = *(pCurrent+0)+k1[0]*c[1]+*(k3+0)*c[3]+*(k4+0)*c[4]+*(k6+0)*c[6];
	}
	barrier.wait();
}
