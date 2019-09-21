#pragma once

#include <cmath>
#include <vector>
#include <complex>

inline double sqr(double x)
{
	return x*x;
}

inline double sign(double x)
{
	if (x<0.0)
		return -1.0;
	return 1.0;
}

inline double complexNorm2Sum(const std::vector<std::complex<double> > &data)
{
	double sum = 0.0;
	for(size_t i=0;i<data.size();++i)
		sum += std::norm(data[i]);
	return sum;
}

inline void meanStdDev(const std::vector<double> &data, double *pMean, double *pStd)
{
	*pMean = 0.0;
	for(const double &x : data)
		*pMean += x;
	*pMean /= data.size();

	*pStd = 0.0;
	for(const double &x : data)
		*pStd += sqr(x-*pMean);
	*pStd = sqrt(*pStd/(data.size()-1));
}

struct LinearRegression
{
	double avgX, stdX;
	double avgY, stdY;
	double a, ua;
	double b, ub;
};

inline void linearRegression(const std::vector<double>& x, const std::vector<double>& y, LinearRegression &lr)
{
	double n = x.size();
	if(x.size()!=y.size())
	{
		// error!
	}

	meanStdDev(x,&lr.avgX,&lr.stdX);
	meanStdDev(y,&lr.avgY,&lr.stdY);

	double numerator = 0.0;
	double denominator = 0.0;

	for(size_t i=0;i<x.size();++i)
	{
		numerator += (x[i] - lr.avgX) * (y[i] - lr.avgY);
		denominator += (x[i] - lr.avgX) * (x[i] - lr.avgX);
	}

	if(fabs(denominator)<1e-16)
	{
		// error!
	}

	lr.a = numerator/denominator;
	lr.b = lr.avgY-lr.a*lr.avgX;

	numerator = 0.0;
	for(size_t i=0;i<x.size();++i)
		numerator += sqr(y[i]-lr.b-lr.a*x[i]);
	lr.ua = sqrt(numerator/denominator/(double)(x.size()-2));

	numerator = 0.0;
	for(size_t i=0;i<x.size();++i)
		numerator += x[i]*x[i];
	lr.ub = lr.ua*sqrt(numerator/(double)x.size());
}
