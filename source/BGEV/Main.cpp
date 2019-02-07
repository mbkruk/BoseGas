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

	std::vector<double> avgs, flucs;

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
			params.batchCount = std::stoul(arg[1]);
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
			std::cerr << "unknown attriute `" << attr << "`" << std::endl;
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

	evolution.create(params);
	evolution.stdinICInit();
	evolution.printParameters();
	for(int i=0;i<params.batchCount;++i)
	{
		evolution.icInit();	
		evolution.evolve(params.batchSize);
		avgs.push_back(evolution.averageNZero());
		flucs.push_back(evolution.fluctuationsNZero(avgs[i]));
		std::cerr  << "Batch " << i+1 << std::setw(dist) << std::setprecision(16) << avgs[i] << std::setw(dist) << 
		flucs[i] << std::setw(dist) << evolution.kineticEnergy()+evolution.potentialEnergy() <<  std::setw(dist) <<
		evolution.nAll() << std::setw(dist) << evolution.momentum() << '\n';
	}
	evolution.destroy();

	return 0;
}
