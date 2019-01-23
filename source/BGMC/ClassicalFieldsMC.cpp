#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <random>
#include <iostream>

#include "ClassicalFieldsMC.hpp"
//#include "Math.hpp"

double ClassicalFieldsMC::momentum()
{
	double m = 0.0;
	for (int32_t n=-nMax;n<=nMax;++n)
		m += n*std::norm(alpha.alpha[n+nMax]);
	return m;
}

void ClassicalFieldsMC::energy(Energy &e)
{
	e = alpha;
}

uint32_t ClassicalFieldsMC::excitedStatesOccupation()
{
	return std::min(static_cast<uint32_t>((N+1.0)/N*std::max(N-std::norm(alpha.alpha[nMax]),0.0)),static_cast<uint32_t>(N+0.5));
}

int32_t ClassicalFieldsMC::steps(std::mt19937 &random, uint32_t count)
{
	std::uniform_real_distribution<double> d1(0.0,1.0);
	std::uniform_real_distribution<double> d2pi(0.0,2.0*M_PI);
	double angle, length;
	int32_t accepted = 0;

	while (count--)
	{
		do
		{
			for (uint32_t i=0;i<alpha0.alpha.size();++i)
			{
				angle = d2pi(random);
				length = delta*d1(random);
				alpha0.alpha[i].real(alpha.alpha[i].real()+length*cos(angle));
				alpha0.alpha[i].imag(alpha.alpha[i].imag()+length*sin(angle));
			}
			alpha0.normalize(N);
		}
		while (fabs(alpha0.magnitude()/N-1.0)>0.01);

		alpha0.update(pInteraction);

		double p = exp(-beta*(alpha0.totalEnergy-alpha.totalEnergy));
		double r = d1(random);

		if (r<p)
		{
			std::swap(alpha,alpha0);
			delta *= 1.0001;
			++accepted;
		}
		else
			delta /= 1.0001;
	}
	return accepted;
}

void ClassicalFieldsMC::generate(std::mt19937 &random)
{
	std::uniform_real_distribution<double> d1(0.0,1.0);
	std::uniform_real_distribution<double> d2pi(0.0,2.0*M_PI);
	std::vector<std::complex<double> > &v = alpha.alpha;
	double length, angle;
	for (int32_t i=0;i<nMax;++i)
	{
		length = d1(random);
		angle = d2pi(random);
		v[i].real(length*cos(angle));
		v[i].imag(length*sin(angle));
	}
	alpha.normalize(N);
	alpha.update(pInteraction);
	delta = 0.05*sqrt(N);
}

void ClassicalFieldsMC::initialize(const BGMCParameters &params)
{
	N = params.particleCount;
	beta = params.beta;
	gamma = params.gamma;
	nMax = params.nMax;

	alpha.alpha.resize(2*nMax+1);
	alpha0.alpha.resize(2*nMax+1);

	pInteraction->prepare(gamma,alpha.alpha.size());

	alpha.asymmetry.resize(nMax);
	alpha0.asymmetry.resize(nMax);
}

void ClassicalFieldsMC::release()
{
}

ClassicalFieldsMC::ClassicalFieldsMC()
{
}

ClassicalFieldsMC::~ClassicalFieldsMC()
{
}

std::vector<std::complex<double> > ClassicalFieldsMC::alphaCopy()
{
	return alpha.alpha;
}

double ClassicalFieldsMC::Alpha::magnitude()
{
	double S = 0.0;
	for (const auto &z : alpha)
		S += std::norm(z);
	return S;
}

void ClassicalFieldsMC::Alpha::normalize(double N)
{
	double sum = 0.0;
	for (const auto &c : alpha)
		sum += std::norm(c);
	double mul = sqrt(N/sum);
	for (auto &c : alpha)
		c *= mul;
}

void ClassicalFieldsMC::Alpha::update(Interaction *pInteraction)
{
	totalEnergy = 0.0;
	kineticEnergy = 0.0;
	interactionEnergy = 0.0;

	int32_t nMax = (static_cast<int32_t>(alpha.size())-1)/2;

	for (int32_t n=-nMax;n<=nMax;++n)
		kineticEnergy += n*n*std::norm(alpha[n+nMax]);

	interactionEnergy = pInteraction->energy(alpha);

	totalEnergy = kineticEnergy+interactionEnergy;

	for (int32_t i=0;i<nMax;++i)
		asymmetry[i] = std::norm(alpha[i])-std::norm(alpha[i+1+nMax]);
}

void ClassicalFieldsMC::NoInteraction::prepare(double gamma_, size_t alphaSize)
{
}

double ClassicalFieldsMC::NoInteraction::energy(const std::vector<std::complex<double> > &alpha) const
{
	return 0.0;
}

void ClassicalFieldsMC::ContactInteraction::prepare(double gamma_, size_t alphaSize)
{
	gamma = gamma_;
	indices.clear();
	for (uint32_t n=0;n<alphaSize;++n)
	for (uint32_t m=0;m<alphaSize;++m)
	for (uint32_t i=0;i<alphaSize;++i)
	for (uint32_t j=0;j<alphaSize;++j)
	if (n+m==i+j)
	{
		indices.push_back({n,m,i,j});
	}
}

double ClassicalFieldsMC::ContactInteraction::energy(const std::vector<std::complex<double> > &alpha) const
{
	double V = 0.0f;
	for (const VI &vi : indices)
	{
		auto c = (alpha[vi.n]*alpha[vi.m]*std::conj(alpha[vi.i])*std::conj(alpha[vi.j]));
		V += std::real(c);
	}
	return gamma*V;
}

void ClassicalFieldsMC::CustomInteraction::prepare(double gamma_, size_t alphaSize)
{
	gamma = gamma_;
	if (coefficients.size()<2*alphaSize-1)
		std::cerr << "Custom Interaction error" << std::endl;
	indices.clear();
	for (uint32_t n=0;n<alphaSize;++n)
	for (uint32_t m=0;m<alphaSize;++m)
	for (uint32_t i=0;i<alphaSize;++i)
	for (uint32_t j=0;j<alphaSize;++j)
	if (n+m==i+j)
	{
		indices.push_back({n,m,i,j});
	}
}

double ClassicalFieldsMC::CustomInteraction::energy(const std::vector<std::complex<double> > &alpha) const
{
	double V = 0.0f;
	for (const VI &vi : indices)
	{
		auto c = (alpha[vi.n]*alpha[vi.m]*std::conj(alpha[vi.i])*std::conj(alpha[vi.j]));
		V += std::real(c)*coefficients[vi.n+coefficients.size()/2-vi.i];
	}
	return gamma*V;
}

ClassicalFieldsMC::CustomInteraction::CustomInteraction(const std::vector<double> &coefficients_)
{
	coefficients.resize(2*coefficients_.size()-1);
	coefficients[coefficients_.size()-1] = coefficients_[0];
	for (size_t i=1;i<coefficients_.size();++i)
	{
		coefficients[coefficients_.size()-1+i] = coefficients_[i];
		coefficients[coefficients_.size()-1-i] = coefficients_[i];
	}
}
