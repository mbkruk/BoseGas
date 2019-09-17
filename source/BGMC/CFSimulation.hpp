#ifndef BGMC_CFSIMULATION_HPP_
#define BGMC_CFSIMULATION_HPP_

#include <cinttypes>
#include <random>
#include <functional>

#include "BGMC.hpp"
#include "ClassicalFieldsMC.hpp"

class CFSimulation
{
public:

	BGMCParameters params;
	std::mt19937 random;
	ClassicalFieldsMC cfmc;
	SimulationInfo batchInfo, totalInfo;
	std::vector<AlphaRecord> alphas;
	BGMC::Energy energy;
	std::vector<double> occupation;
	
	void batch(bool collect);

	void finish();
	
	static ClassicalFieldsMC::Interaction* createInteraction(const BGMCParameters &params_);

	int32_t initialize(const BGMCParameters &params_, std::function<void(CFSimulation&sim)> initCallback, std::function<void(CFSimulation&sim)> acceptCallback);
	void release();
};

int32_t bgSingleCFSimulation(BGMCParameters &params);
int32_t bgCFCompare(const BGMCParameters &params);

#endif
