#pragma once

#include <cinttypes>
#include <complex>

#include "BGMC.hpp"

class ClassicalFieldsMC : public BGMC
{
public:

	class Interaction
	{
	private:
	public:
		virtual void prepare(double gamma_, size_t alphaSize) = 0;
		virtual double energy(const std::vector<std::complex<double> > &alpha) const = 0;
	};

	class NoInteraction : public Interaction
	{
	public:
		void prepare(double gamma_, size_t alphaSize) override;
		double energy(const std::vector<std::complex<double> > &alpha) const override;
	};

	class ContactInteraction : public Interaction
	{
	private:
		double gamma;
		struct VI
		{
			uint32_t n;
			uint32_t m;
			uint32_t i;
			uint32_t j;
		};
		std::vector<VI> indices;
	public:
		void prepare(double gamma_, size_t alphaSize) override;
		double energy(const std::vector<std::complex<double> > &alpha) const override;
	};

	class CustomInteraction : public Interaction
	{
	private:
		double gamma;
		std::vector<double> coefficients;
		struct VI
		{
			uint32_t n;
			uint32_t m;
			uint32_t i;
			uint32_t j;
		};
		std::vector<VI> indices;
	public:
		void prepare(double gamma_, size_t alphaSize) override;
		double energy(const std::vector<std::complex<double> > &alpha) const override;
		CustomInteraction(const std::vector<double> &coefficients_);
	};

private:

	double N;
	double beta;
	double gamma;

	double deltaCoeff;
	double delta;
	double deltaDelta;
	bool constDelta;

	int32_t baseNMax;
	int32_t nMax;

	struct Alpha : public Energy
	{
		std::vector<double> asymmetry;
		std::vector<std::complex<double> > alpha;

		double magnitude();
		void normalize(double N);
		void update(Interaction *pInteraction);
	};

	Alpha alpha, alpha0;

	Interaction *pInteraction = nullptr;

public:

	static double getOptimalContactGamma(BGMCParameters &params);

	inline double getInitialDetla() const
	{
		return deltaCoeff*sqrt(2*baseNMax+1)/(double)(1*(nMax-baseNMax)+1)/sqrt(1.0+std::abs(N*gamma));
	}

	inline double getDelta() const
	{
		return delta;
	}
	
	inline void setDelta(double d)
	{
		delta = d;
	}

	inline double getRelativeDelta() const
	{
		return delta/getInitialDetla();
	}

	inline const Alpha& currentAlpha()
	{
		return alpha;
	}

	std::vector<std::complex<double> > alphaCopy();

	inline int32_t getNMax() const
	{
		return nMax;
	}

	/*inline int32_t setNMax(int32_t newNMax)
	{
		nMa = newNMax;
	}*/

	inline void setInteraction(Interaction *pNewInteraction)
	{
		pInteraction = pNewInteraction;
	}
	
	inline Interaction* getInteraction()
	{
		return pInteraction;
	}

	double momentum() override;
	void energy(Energy &e) override;

	uint32_t excitedStatesOccupation() override;
	double groundStateOccupation();

	void getOccupation(std::vector<double> &occupation);

	int32_t steps(std::mt19937 &random, uint32_t count) override;

	void generate(std::mt19937 &random) override;

	void initialize(const BGMCParameters &params) override;
	void release() override;

	ClassicalFieldsMC();
	virtual ~ClassicalFieldsMC();
};
