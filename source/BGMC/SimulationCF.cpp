#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>

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
	double pAcceptMean, pAcceptStdDev, pAcceptMeanStdDev;
	std::vector<double> n0;
	double n0Mean, n0StdDev, n0MeanStdDev;
	std::vector<double> asymmetry;
	double asymmetryMean, asymmetryStdDev;
	std::vector<double> energy;
	double energyMean, energyStdDev, energyMeanStdDev;
	std::vector<double> momentum;
	double momentumMean, momentumStdDev, momentumMeanStdDev;

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
		pAcceptMeanStdDev = pAcceptStdDev/sqrt(pAccept.size()-1);
		meanStdDev(n0,&n0Mean,&n0StdDev);
		n0MeanStdDev = n0StdDev/sqrt(n0.size()-1);
		meanStdDev(asymmetry,&asymmetryMean,&asymmetryStdDev);
		meanStdDev(energy,&energyMean,&energyStdDev);
		energyMeanStdDev = energyStdDev/sqrt(energy.size()-1);
		meanStdDev(momentum,&momentumMean,&momentumStdDev);
		momentumMeanStdDev = momentumStdDev/sqrt(momentum.size()-1);
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
			<< " " << std::setw(15) << "n0 mean" << std::setw(15) << "+/-" << std::setw(15) << "n0 std dev"
			//<< " " << std::setw(15) << "asym mean" << std::setw(15) << "asym std dev"
			<< " " << std::setw(15) << "E mean" << std::setw(15) << "+/-" << std::setw(15) << "E std dev"
			<< " " << std::setw(15) << "P mean" << std::setw(15) << "+/-" << std::setw(15) << "P std dev"
			<< std::endl;
	}

	template<typename X> void print(std::ostream &os, X batch)
	{
		os << std::setw(8) << batch
			<< " " << std::setw(15) << pAcceptMean
			<< " " << std::setw(15) << n0Mean << std::setw(15) << n0MeanStdDev << std::setw(15) << n0StdDev
			//<< " " << std::setw(15) << asymmetryMean << std::setw(15) << asymmetryStdDev
			<< " " << std::setw(15) << energyMean << std::setw(15) << energyMeanStdDev << std::setw(15) << energyStdDev
			<< " " << std::setw(15) << momentumMean << std::setw(15) << momentumMeanStdDev << std::setw(15) << momentumStdDev
			<< std::endl;
	}
};

int32_t bgSimulationCF(BGMCParameters &params)
{
	std::mt19937 random;
	ClassicalFieldsMC cfmc;
	ClassicalFieldsMC::Interaction *pInteraction = nullptr;
	Info batchInfo, totalInfo;
	std::vector<AlphaRecord> alphas;

	if (params.gamma!=0.0 && params.interactionType!="none")
	{
		if (params.interactionType=="contact")
		{
			//std::cerr << "using contact interaction" << std::endl;
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
	{
		pInteraction = new ClassicalFieldsMC::NoInteraction;
		//std::cerr << "no interaction" << std::endl;
	}
	cfmc.setInteraction(pInteraction);
	random.seed(params.seed);

	cfmc.initialize(params);

	Info::printHead(std::cout);

	BGMC::Energy energy;

	totalInfo.clear();

	do
	{
		cfmc.generate(random);
		batchInfo.clear();
		for (uint32_t i=0;i<params.batchSize;++i)
		{
			double a = cfmc.steps(random,1);
			cfmc.energy(energy);
			double n0 = cfmc.groundStateOccupation();
			double p = cfmc.momentum();
			batchInfo.append(a,n0,0.0,energy.totalEnergy,p);
		}
		batchInfo.process();
		batchInfo.print(std::cout,"init");
	}
	while (std::abs(batchInfo.momentumMean)/batchInfo.momentumMeanStdDev>=3.0);

	while (
		std::abs(batchInfo.pAcceptMean-0.5)/batchInfo.pAcceptMeanStdDev>=4.0
		|| std::abs(batchInfo.momentumMean)/batchInfo.momentumMeanStdDev>=3.0
		)
	{
		batchInfo.clear();
		for (uint32_t i=0;i<params.batchSize;++i)
		{
			double a = cfmc.steps(random,1);
			cfmc.energy(energy);
			double n0 = cfmc.groundStateOccupation();
			double p = cfmc.momentum();
			batchInfo.append(a,n0,0.0,energy.totalEnergy,p);
		}
		batchInfo.process();
		batchInfo.print(std::cout,"accept");
	}

	for (uint32_t batch=1;batch<params.batchCount;++batch)
	{
		batchInfo.clear();
		for (uint32_t i=0;i<params.batchSize;++i)
		{
			double a = cfmc.steps(random,1);
			cfmc.energy(energy);
			double n0 = cfmc.groundStateOccupation();
			double p = cfmc.momentum();
			batchInfo.append(a,n0,0.0,energy.totalEnergy,p);
			if ((i+1)%params.skip==0)
			{
				totalInfo.append(a,n0,0.0,energy.totalEnergy,p);
				alphas.push_back({energy.totalEnergy,p,0.0,0.0,0.0,cfmc.alphaCopy()});
			}
		}
		batchInfo.process();
		batchInfo.print(std::cout,batch);
	}

	totalInfo.process();
	totalInfo.print(std::cout,"total");

	for (auto &ar : alphas)
	{
		ar.Ed = (ar.E-totalInfo.energyMean)/totalInfo.energyStdDev;
		ar.Pd = ar.P/totalInfo.momentumStdDev;
		ar.distance = sqrt(sqr(ar.Ed)+sqr(ar.Pd));
	}

	std::sort(alphas.begin(),alphas.end(),[](const AlphaRecord &a, const AlphaRecord &b){return a.distance<b.distance;});

	std::fstream output_file;
	output_file.open(params.output,std::ios::out);

	output_file << params.particleCount << '\n';
	output_file << params.nMax+params.extraModePairs << '\n';
	output_file << params.gamma << '\n';
	output_file << params.interactionType << '\n' << '\n';

	if (params.interactionType=="custom")
	{
		for (int_fast32_t i=0;i<2*(params.nMax+params.extraModePairs)+1;++i)
			output_file << params.interactionCoefficients[i] << '\n';
		output_file << '\n';
	}

	for (int_fast32_t i=0; i<2*(params.nMax+params.extraModePairs)+1;++i)
		output_file << std::real(alphas[0].alpha[i]) << '\n';

	for (int_fast32_t i=0; i<2*(params.nMax+params.extraModePairs)+1;++i)
		output_file << std::imag(alphas[0].alpha[i]) << '\n';

	output_file << '\n';
	Info::printHead(output_file);
	totalInfo.print(output_file,"total");

	output_file.close();

	cfmc.release();
	delete pInteraction;
	return 0;
}
