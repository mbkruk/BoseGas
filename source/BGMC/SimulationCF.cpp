#include <iostream>
#include <iomanip>

#include "BGMC.hpp"
#include "ClassicalFieldsMC.hpp"
#include "../BGCommon/Math.hpp"

struct AlphaRecord
{
	double E;
	double P;
	double Ed, Pd;
	double distance;
	std::vector<std::complex<double> > alpha;
};

struct Info
{
	std::vector<double> pAccept;
	double pAcceptMean, pAcceptStdDev;
	std::vector<double> n0;
	double n0Mean, n0StdDev;
	std::vector<double> asymmetry;
	double asymmetryMean, asymmetryStdDev;
	std::vector<double> energy;
	double energyMean, energyStdDev;
	std::vector<double> momentum;
	double momentumMean, momentumStdDev;

	void append(double accepted, double N0, double A, double E, double P)
	{
		pAccept.push_back(accepted);
		n0.push_back(N0);
		asymmetry.push_back(A);
		energy.push_back(E);
		momentum.push_back(P);
	}

	void process()
	{
		meanStdDev(pAccept,&pAcceptMean,&pAcceptStdDev);
		meanStdDev(n0,&n0Mean,&n0StdDev);
		meanStdDev(asymmetry,&asymmetryMean,&asymmetryStdDev);
		meanStdDev(energy,&energyMean,&energyStdDev);
		meanStdDev(momentum,&momentumMean,&momentumStdDev);
	}

	void clear()
	{
		pAccept.clear();
		n0.clear();
		asymmetry.clear();
		energy.clear();
		momentum.clear();
	}

	static void printHead(std::ostream &os)
	{
		os << std::setw(8) << "batch"
			<< " " << std::setw(15) << "P accept"
			<< " " << std::setw(15) << "n0 mean" << std::setw(15) << "n0 std dev"
			<< " " << std::setw(15) << "asym mean" << std::setw(15) << "asym std dev"
			<< " " << std::setw(15) << "E mean" << std::setw(15) << "E std dev"
			<< " " << std::setw(15) << "P mean" << std::setw(15) << "P std dev"
			<< std::endl;
	}

	void print(std::ostream &os, uint32_t batchIndex)
	{
		os << std::setw(8) << batchIndex
			<< " " << std::setw(15) << pAcceptMean
			<< " " << std::setw(15) << n0Mean << std::setw(15) << n0StdDev
			<< " " << std::setw(15) << asymmetryMean << std::setw(15) << asymmetryStdDev
			<< " " << std::setw(15) << energyMean << std::setw(15) << energyStdDev
			<< " " << std::setw(15) << momentumMean << std::setw(15) << momentumStdDev
			<< std::endl;
	}
};

int32_t bgSimulationCF(const BGMCParameters &params)
{
	std::mt19937 random;
	ClassicalFieldsMC cfmc;
	ClassicalFieldsMC::Interaction *pInteraction = nullptr;
	Info batchInfo;

	if (params.gamma!=0.0 && params.interactionType!="none")
	{
		if (params.interactionType=="contact")
		{
			std::cerr << "using contact interaction" << std::endl;
			pInteraction = new ClassicalFieldsMC::ContactInteraction;
		}
		else
		if (params.interactionType=="custom")
		{
			pInteraction = new ClassicalFieldsMC::CustomInteraction(params.interactionCoefficients);
		}
		else
		{
			std::cerr << "unknown interaction type `" << params.interactionType << "`" << std::endl;
			return 1;
		}
	}
	else
		pInteraction = new ClassicalFieldsMC::NoInteraction;
	cfmc.setInteraction(pInteraction);
	random.seed(params.seed);

	cfmc.initialize(params);
	cfmc.generate(random);

	Info::printHead(std::cerr);

	BGMC::Energy energy;

	for(uint32_t batch=0;batch<params.batchCount;++batch)
	{
		for (uint32_t i=0;i<params.batchSize;++i)
		{
			double p = cfmc.steps(random,1);
			cfmc.energy(energy);
			batchInfo.append(p,cfmc.groundStateOccupation(),0.0,energy.totalEnergy,cfmc.momentum());
		}
		batchInfo.process();
		batchInfo.print(std::cerr,batch);
	}

	cfmc.release();
	delete pInteraction;
	return 0;
}
