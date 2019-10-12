#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <thread>
#include <atomic>
#include <mutex>

#include "CFSimulation.hpp"
#include "../BGCommon/Math.hpp"
#include "../BGCommon/Barrier.hpp"
#include "../BGCommon/Random.hpp"

ClassicalFieldsMC::Interaction* CFSimulation::createInteraction(const BGMCParameters &params)
{
	ClassicalFieldsMC::Interaction *pInteraction = nullptr;

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
			return nullptr;
		}
	}
	else
	{
		pInteraction = new ClassicalFieldsMC::NoInteraction;
		//std::cerr << "no interaction" << std::endl;
	}

	return pInteraction;
}

inline double regressionScore(const BGMCParameters &params, const LinearRegression &lr, const SimulationInfo &batchInfo)
{
	return (lr.a*params.batchSize+lr.b-batchInfo.energyMean)/sqrt(sqr(lr.ua*params.batchSize)+sqr(lr.ub)+sqr(batchInfo.energyMeanStdDev));
}

static int32_t optSkip(const BGMCParameters &params, const CFSimulation &sim)
{
	double D = (2.0/M_PI)*sim.deltaMean*sqrt(2*(params.nMax+params.extraModePairs)+1);
	double N = M_PI*sqrt((double)params.particleCount);
	return std::max((int32_t)(N/std::max(D,N/1024.0)+0.5),4);
}

int32_t CFSimulation::initialize(const BGMCParameters &params_, std::function<void(CFSimulation&sim)> initCallback, std::function<void(CFSimulation&sim)> prepareCallback)
{
	params = params_;

	random.seed(params.seed);

	cfmc.initialize(params);

	totalInfo.clear(cfmc.getNMax());

	occupation.resize(2*cfmc.getNMax()+1);

	batchRound.clear();
	for (uint32_t i=0;i<params.batchSize;++i)
		batchRound.push_back((double)i-0.5*(params.batchSize-1));

	double lrscore;

	uint32_t cnt = 0;
	do
	{
		delta.clear();
		if ((cnt%2)==0)
			cfmc.generate(random);
		batchInfo.clear(cfmc.getNMax());
		for (uint32_t i=0;i<params.batchSize*params.skip;++i)
		if (((i+1)%params.skip)==0)
		{
			delta.push_back(cfmc.getDelta());
			double a = cfmc.steps(random,1);
			cfmc.energy(energy);
			cfmc.getOccupation(occupation);
			double p = cfmc.momentum();
			batchInfo.append(a,occupation,0.0,energy.totalEnergy,p);
		}
		else
			cfmc.steps(random,1);
		meanStdDev(delta,&deltaMean,&deltaStdDev);
		deltaMeanStdDev = deltaStdDev/sqrt(delta.size()-1);
		batchInfo.process();
		linearRegression(batchRound,batchInfo.energy,lr);
		lrscore = regressionScore(params,lr,batchInfo);
		initCallback(*this);
		++cnt;
	}
	while (!params.acceptTest
		&& (
			(
				params.useConstDelta && batchInfo.pAcceptMeanStdDev>0.01)
				|| (!params.useConstDelta && std::abs(batchInfo.pAcceptMean-0.5)/batchInfo.pAcceptMeanStdDev>=3.0)
			)
		);

	params.skip = optSkip(params,*this);
	std::cerr << "set skip " << params.skip << std::endl;

	while (
		(params.useConstDelta && batchInfo.pAcceptMeanStdDev>0.01)
		|| (!params.useConstDelta && std::abs(batchInfo.pAcceptMean-0.5)/batchInfo.pAcceptMeanStdDev>=3.0)
		|| std::abs(batchInfo.momentumMean)/batchInfo.momentumMeanStdDev>=3.0
		|| fabs(lrscore)>=3.0
		)
	{
		delta.clear();
		batchInfo.clear(cfmc.getNMax());
		for (uint32_t i=0;i<params.batchSize*params.skip;++i)
		if (((i+1)%params.skip)==0)
		{
			delta.push_back(cfmc.getDelta());
			double a = cfmc.steps(random,1);
			cfmc.energy(energy);
			cfmc.getOccupation(occupation);
			double p = cfmc.momentum();
			batchInfo.append(a,occupation,0.0,energy.totalEnergy,p);
		}
		else
			cfmc.steps(random,1);
		meanStdDev(delta,&deltaMean,&deltaStdDev);
		deltaMeanStdDev = deltaStdDev/sqrt(delta.size()-1);
		batchInfo.process();
		linearRegression(batchRound,batchInfo.energy,lr);
		lrscore = regressionScore(params,lr,batchInfo);
		prepareCallback(*this);
	}

	return 0;
}

void CFSimulation::batch(bool collect)
{
	batchInfo.clear(cfmc.getNMax());
	for (uint32_t i=0;i<params.batchSize*params.skip;++i)
	{
		double a = cfmc.steps(random,1);
		cfmc.energy(energy);
		double p = cfmc.momentum();
		if ((i+1)%params.skip==0)
		{
			cfmc.getOccupation(occupation);
			batchInfo.append(a,occupation,0.0,energy.totalEnergy,p);
			if (collect)
			{
				totalInfo.append(a,occupation,0.0,energy.totalEnergy,p);
				if (params.skipOutput && (i/params.skip+1)%params.skipOutput==0)
					alphas.push_back({energy.totalEnergy,p,0.0,0.0,0.0,cfmc.alphaCopy()});
			}
		}
	}
	batchInfo.process();
}

void CFSimulation::finish()
{
	totalInfo.process();

	for (auto &ar : alphas)
	{
		ar.Ed = (ar.E-totalInfo.energyMean)/totalInfo.energyStdDev;
		ar.Pd = ar.P/totalInfo.momentumStdDev;
		ar.distance = sqrt(sqr(ar.Ed)+sqr(ar.Pd));
	}

	if (params.sort)
		std::sort(alphas.begin(),alphas.end(),[](const AlphaRecord &a, const AlphaRecord &b){return a.distance<b.distance;});
}

void CFSimulation::release()
{
	cfmc.release();
}

int32_t bgSingleCFSimulation(BGMCParameters &params)
{
	CFSimulation sim;

	SimulationInfo::printHead(std::cout);

	sim.cfmc.setInteraction(CFSimulation::createInteraction(params));
	sim.initialize(params,
		[](CFSimulation &sim){sim.batchInfo.print(std::cout,"init");},
		[](CFSimulation &sim){sim.batchInfo.print(std::cout,"prepare");}
	);

	if (params.acceptTest)
	{
		std::cout << "delta = " << sim.deltaMean << " +/- " << sim.deltaMeanStdDev << std::endl;
		std::cout << "optimal skip parameter: " << optSkip(params,sim) << std::endl;
		return 0;
	}

	for (uint32_t batch=1;batch<=params.batchCount;++batch)
	{
		sim.batch(batch>=params.firstBatch);
		sim.batchInfo.print(std::cout,batch);
	}

	sim.finish();

	sim.totalInfo.print(std::cout,"total");
	std::cout << std::endl;
	sim.totalInfo.printOccupation(std::cout);

	if (!params.output.empty())
	{
		std::fstream output_file;
		output_file.open(params.output,std::ios::out);

		output_file << params.particleCount << '\n';
		output_file << params.nMax+params.extraModePairs << '\n';
		output_file << params.gamma << '\n';
		output_file << params.interactionType << '\n';

		if (params.outputStyle=="modified")
			output_file << sim.alphas.size() << '\n' << '\n';
		else
			output_file << "1" << '\n' << '\n';

		if (params.outputStyle=="modified")
		{
			for (int_fast32_t i=0; i<sim.alphas.size();++i)
			{
				for (int_fast32_t j=0; j<2*(params.nMax+params.extraModePairs)+1;++j)
					output_file << std::real(sim.alphas[i].alpha[j]) << " ";

				for (int_fast32_t j=0; j<2*(params.nMax+params.extraModePairs)+1;++j)
					output_file << std::imag(sim.alphas[i].alpha[j]) << " ";

				output_file << '\n';
			}

			output_file << '\n';
		}
		else
		{
			if (params.interactionType=="custom")
			{
				for (int_fast32_t i=0;i<2*(params.nMax+params.extraModePairs)+1;++i)
					output_file << params.interactionCoefficients[i] << '\n';
				output_file << '\n';
			}

			for (int_fast32_t i=0; i<2*(params.nMax+params.extraModePairs)+1;++i)
				output_file << std::real(sim.alphas[0].alpha[i]) << '\n';

			for (int_fast32_t i=0; i<2*(params.nMax+params.extraModePairs)+1;++i)
				output_file << std::imag(sim.alphas[0].alpha[i]) << '\n';

			output_file << '\n';
		}

		SimulationInfo::printHead(output_file);
		sim.totalInfo.print(output_file,"total");
		output_file << std::endl;
		sim.totalInfo.printOccupation(output_file);

		output_file.close();
	}

	sim.release();
	delete sim.cfmc.getInteraction();
	return 0;
}

void bgCFCompareThread(uint32_t index, const BGMCParameters &params_, CFSimulation *pSims, Barrier &barrier, std::atomic_bool &run, std::mutex &logMutex)
{
	CFSimulation &sim = pSims[index];
	BGMCParameters params = params_;
	if (index!=0)
		params.seed = generateSeed();
	params.extraModePairs += (int32_t)index-1;

	sim.cfmc.setInteraction(CFSimulation::createInteraction(params));
	sim.initialize(params,
		[&](CFSimulation &sim2)
		{
			std::lock_guard<std::mutex> lk(logMutex);
			std::cerr << std::setw(8) << index;
			sim.batchInfo.print(std::cerr,"init");
		},
		[&](CFSimulation &sim2)
		{
			std::lock_guard<std::mutex> lk(logMutex);
			std::cerr << std::setw(8) << index;
			sim.batchInfo.print(std::cerr,"prepare");
		}
	);

	barrier.wait();

	for (uint32_t batch=1;run.load() && batch<=params.batchCount;++batch)
	{
		sim.batch(batch>=params.firstBatch);
		sim.totalInfo.process();

		barrier.wait();

		if (index==0)
		{
			if (!sim.totalInfo.energy.empty())
			{
				for (uint32_t i=0;i<3;++i)
				{
					std::lock_guard<std::mutex> lk(logMutex);
					std::cerr << std::setw(7) << i;
					pSims[i].totalInfo.print(std::cerr,batch-params.firstBatch+1);
				}

				double ldiff = pSims[1].totalInfo.energyMean-pSims[0].totalInfo.energyMean;
				double uldiff = sqrt(sqr(pSims[1].totalInfo.energyMeanStdDev)+sqr(pSims[0].totalInfo.energyMeanStdDev)+1e-24);
				double rdiff = pSims[2].totalInfo.energyMean-pSims[1].totalInfo.energyMean;
				double urdiff = sqrt(sqr(pSims[2].totalInfo.energyMeanStdDev)+sqr(pSims[1].totalInfo.energyMeanStdDev)+1e-24);
				double l = ldiff/uldiff;
				double r = rdiff/urdiff;
				std::cout << l << " " << r << std::endl;
				if (fabs(l)>3.0 && fabs(r)>3.0)
					run.store(false);
			}
		}

		barrier.wait();
	}
	delete sim.cfmc.getInteraction();
}

int32_t bgCFCompare(const BGMCParameters &params)
{
	std::atomic_bool run(true);
	Barrier barrier(3);
	std::mutex logMutex;

	CFSimulation sim[3];

	std::cerr << std::setw(8) << "thread";
	SimulationInfo::printHead(std::cerr);

	std::thread th[2];
	th[1] = std::thread(bgCFCompareThread,2,std::cref(params),sim,std::ref(barrier),std::ref(run),std::ref(logMutex));
	th[0] = std::thread(bgCFCompareThread,1,std::cref(params),sim,std::ref(barrier),std::ref(run),std::ref(logMutex));

	bgCFCompareThread(0,params,sim,barrier,run,logMutex);

	th[0].join();
	th[1].join();
	return 0;
}
