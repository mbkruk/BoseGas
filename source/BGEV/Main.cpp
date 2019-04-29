#include <iostream>
#include <cinttypes>
#include <iomanip>
#include <vector>

#include "../BGCommon/ProgramOptions.hpp"
#include "../BGCommon/ConfigFile.hpp"
#include "BGEV.hpp"

const int_fast32_t dist = 2*sizeof(double)+9;

int main(int argc, const char *argv[])
{
	ProgramOptions po;
	ConfigFile cf;
	int r;
	bool help = false;
	std::vector<std::string> inputFiles;
	BGEVParameters params;
	BGEvolution evolution;

	params.particleCount = 1;
	params.nMax = 1;
	params.gamma = 1;
	params.batchSize = 10000;
	params.batchCount = 1;
	params.threadCount = 0;

	po.addOption("-h","--help","produce help message",1,[&](const char *[]){help=true;return 0;});

	po.addOption("-o","--output <file>","set output file",2,
		[&](const char *arg[])
		{
			params.output = arg[1];
			return 0;
		}
	);

	po.addOption("-H","--step-size <number>","set H(step size)",2,
		[&](const char *arg[])
		{
			params.h = std::stod(arg[1]);
			return 0;
		}
	);

	po.addOption("-BS","--Batch-size <count>","set batch size",2,
		[&](const char *arg[])
		{
			params.batchSize = std::stoul(arg[1]);
			return 0;
		}
	);

	po.addOption("-BC","--Batch-count <count>","set batch count",2,
		[&](const char *arg[])
		{
			params.batchCount = std::stoul(arg[1]);
			return 0;
		}
	);

	po.addOption("-t","--thread-count <count>","set extra thread count",2,
		[&](const char *arg[])
		{
			params.threadCount = std::stoul(arg[1]);
			return 0;
		}
	);

	po.addOption("","",nullptr,1,
		[&](const char *arg[])
		{
			if (**arg=='-')
			{
				std::cerr << "unknown option `" << *arg << "`" << std::endl;
				return 2;
			}
			inputFiles.push_back(*arg);
			return 0;
		}
	);

	if (r=po.parseOptions(argc,argv))
		return r;

	if (help)
	{
		std::cerr << "usage: bgev <options> ..." << std::endl;
		std::cerr << "options:" << std::endl;
		po.printDescriptions(std::cerr);
		return 0;
	}

	cf.addAttribute("H",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->h))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("BS",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->batchSize))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("BC",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->batchCount))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("",
		[](std::istream &is, const char *attr, void *pData)
		{
			std::cerr << "unknown attribute `" << attr << "`" << std::endl;
			return 2;
		}
	);

	for (std::string &f : inputFiles)
	{
		if (r=cf.parse(f.c_str(),&params))
			return r;
	}

	std::cin >> params.particleCount;
	std::cin >> params.nMax;
	std::cin >> params.gamma;
	std::cin >> params.interactionType;
	std::cin >> params.alphaCount;

	if (params.interactionType=="custom")
	{
		double x;
		for (int_fast32_t i=0;i<=2*params.nMax;++i)
		{
			std::cin >> x;
			params.reducedCoefficients.push_back(x);
		}
	}
	else
	if (params.interactionType!="contact")
	{
		std::cout << "Unknown interaction type";
		return 1;
	}

	evolution.create(params);
	evolution.stdinInit();
	evolution.printParameters();

	for (int_fast32_t i=0;i<params.batchCount;++i)
	{
		evolution.icInit();
		if (params.threadCount>0)
			evolution.evolve2();
		else
			evolution.evolve();
		evolution.calcAverages();
		std::cerr  << std::setw(dist) << i+1 << std::setw(dist) << std::setprecision(16) << evolution.avgs[i] << std::setw(dist)
			<< evolution.flucs[i] << std::setw(dist) << evolution.kineticEnergy()+evolution.potentialEnergy() <<  std::setw(dist)
			<< evolution.nAll() << std::setw(dist) << evolution.momentum() << '\n';
		evolution.saveToFile(i);
		if (i==params.batchCount-1)
			evolution.lastBatch();
	}
	evolution.destroy();

	return 0;
}
