#ifndef BGMC_HPP_
#define BGMC_HPP_

#include <cinttypes>
#include <string>
#include <random>

static constexpr double kCutoffConst = 0.58;

struct BGMCParameters
{
	uint32_t particleCount;
	int32_t nMax;
	double beta;
	double gamma;
	uint32_t seed;
	uint32_t batchCount, batchSize;
	std::string interactionType;
	std::vector<double> interactionCoefficients;
	std::string output;
};

class BGMC
{
public:

	struct Energy
	{
		double totalEnergy;
		double kineticEnergy;
		double interactionEnergy;
	};

	virtual double momentum() = 0;
	virtual void energy(Energy &e) = 0;

	virtual uint32_t excitedStatesOccupation() = 0;

	virtual int32_t steps(std::mt19937 &random, uint32_t count) = 0;

	virtual void generate(std::mt19937 &random) = 0;

	virtual void initialize(const BGMCParameters &params) = 0;
	virtual void release() = 0;

	virtual ~BGMC();
};

int32_t bgSimulationCF(const BGMCParameters &params);

#endif
