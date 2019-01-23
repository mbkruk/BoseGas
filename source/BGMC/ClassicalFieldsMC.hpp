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

	double delta;

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

	inline double getDelta() const
	{
		return delta;
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

	double momentum() override;
	void energy(Energy &e) override;

	uint32_t excitedStatesOccupation() override;
	double groundStateOccupation();

	int32_t steps(std::mt19937 &random, uint32_t count) override;

	void generate(std::mt19937 &random) override;

	void initialize(const BGMCParameters &params) override;
	void release() override;

	ClassicalFieldsMC();
	virtual ~ClassicalFieldsMC();
};
