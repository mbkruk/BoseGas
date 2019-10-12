#include <iostream>
#include <cinttypes>
#include <string>
#include <vector>
#include <random>
#include <cmath>

#include "../BGCommon/ProgramOptions.hpp"
#include "../BGCommon/ConfigFile.hpp"
#include "../BGCommon/Random.hpp"
#include "BGMC.hpp"
#include "ClassicalFieldsMC.hpp"
#include "CFSimulation.hpp"

int main(int argc, const char *argv[])
{
	ProgramOptions po;
	ConfigFile cf;
	int r;
	bool help = false;
	std::vector<std::string> inputFiles;
	BGMCParameters params;

	params.particleCount = 1;
	params.nMax = 0;
	params.extraModePairs = 0;
	params.beta = 1.0;
	params.gamma = 0.0;
	params.seed = generateSeed();
	params.batchCount = 8;
	params.batchSize = 1024*16;
	params.skip = 4;
	params.skipOutput = 1;
	params.firstBatch = 1;
	params.interactionType = "contact";
	params.useConstDelta = false;
	params.sort = true;
	params.acceptTest = false;
	params.delta = 1.0;
	params.betaRatio = 1.0;

	bool gammaSet = false;
	bool compare = false;

	po.addOption("-h","--help","produce help message",1,[&](const char *[]){help=true;return 0;});

	po.addOption("-o","--output <file>","set output file",2,
		[&](const char *arg[])
		{
			params.output = arg[1];
			return 0;
		}
	);

	po.addOption("-N","--particle-count <count>","set particle count",2,
		[&](const char *arg[])
		{
			params.particleCount = std::stoul(arg[1]);
			return 0;
		}
	);

	po.addOption("-n","--nmax <count>","set Nmax and beta optimized for groud state occupation",2,
		[&](const char *arg[])
		{
			params.nMax = std::stol(arg[1]);
			params.beta = kCutoffConst/params.nMax/params.nMax;
			return 0;
		}
	);

	po.addOption("-w","--hwhm-cutoff <count>","set Nmax>=2 and beta optimized for HWHM",2,
		[&](const char *arg[])
		{
			double mult[] = {
				0.0,
				0.0,
				0.42522181368612244,
				0.6184360448913647,
				0.7626433192054166,
				0.8812398719072574,
				0.9827619383397457,
				1.0708576691904812,
				1.1478793269503433,
				1.215539164651759,
				1.2750995561728224,
				1.3278649752647078,
				1.3753891615680687,
				1.4189174944290668,
				1.4587286304414895,
				1.4940182544970568
			};
			params.nMax = std::stol(arg[1]);
			if (params.nMax<2)
			{
				std::cerr << "invalid Nmax" << std::endl;
				return 1;
			}
			params.betaRatio = mult[params.nMax];
			params.beta = kCutoffConst*mult[params.nMax]/params.nMax/params.nMax;
			return 0;
		}
	);

	po.addOption("-e","--extra <count>","extra mode pairs",2,
		[&](const char *arg[])
		{
			params.extraModePairs = std::stol(arg[1]);
			if (!gammaSet)
				params.gamma = ClassicalFieldsMC::getOptimalContactGamma(params);
			return 0;
		}
	);

	po.addOption("-b","--beta <number>","set beta (inverse temperature)",2,
		[&](const char *arg[])
		{
			params.beta = std::stod(arg[1]);
			params.nMax = static_cast<int32_t>(sqrt(kCutoffConst/params.beta));
			return 0;
		}
	);

	po.addOption("-T","--temperature <number>","set temperature",2,
		[&](const char *arg[])
		{
			params.beta = 1.0/std::stod(arg[1]);
			params.nMax = static_cast<int32_t>(sqrt(kCutoffConst/params.beta));
			return 0;
		}
	);

	po.addOption("-S","--batch-size <count>","set batch size (default is 65536)",2,
		[&](const char *arg[])
		{
			params.batchSize = std::stol(arg[1]);
			return 0;
		}
	);

	po.addOption("-B","--batch-count <count>","set batch count (default is 8)",2,
		[&](const char *arg[])
		{
			params.batchCount = std::stol(arg[1]);
			return 0;
		}
	);

	po.addOption("-g","--gamma <number>","interaction strength",2,
		[&](const char *arg[])
		{
			params.gamma = std::stod(arg[1]);
			gammaSet = true;
			return 0;
		}
	);

	po.addOption("-I","--interaction <file>","set custom interaction",2,
		[&](const char *arg[])
		{
			params.interactionType = "custom";
			std::ifstream f(arg[1]);
			if (!f.is_open())
			{
				std::cerr << "file `" << arg[1] << "` not open" << std::endl;
				return 1;
			}
			{
				std::string word;
				if (!(f>>word) || word!="interaction_coefficients")
				{
					std::cerr << "wrong file format" << std::endl;
					return 1;
				}
			}
			double c;
			while (!f.eof())
			{
				if (!(f>>c))
					break;
				params.interactionCoefficients.push_back(c);
			}

			f.close();
			return 0;
		}
	);

	po.addOption("-a","--alpha-output","output unsorted sets of alphas (for random walk and local fluctuations)",1,
		[&](const char *[])
		{
			params.outputStyle = "modified";
			params.sort = false;
			return 0;
		}
	);

	po.addOption("-s","--skip <s>","use every s-th alpha set to calculate averages and fluctuations (default is 4)",2,
		[&](const char *arg[])
		{
			params.skip = std::stol(arg[1]);
			return 0;
		}
	);

	po.addOption("-O","--skip-output <s>","save every s-th set of alphas (default is 4)",2,
		[&](const char *arg[])
		{
			params.skipOutput = std::stol(arg[1]);
			return 0;
		}
	);

	po.addOption("-d","--delta <number>","set constant MC delta",2,
		[&](const char *arg[])
		{
			params.useConstDelta = true;
			params.delta = std::stod(arg[1]);
			return 0;
		}
	);

	po.addOption("-f","--first-batch <n>","collect alphas only from batch n and later",2,
		[&](const char *arg[])
		{
			params.firstBatch = std::stoi(arg[1]);
			return 0;
		}
	);

	po.addOption("-c","--compare","compare energies",1,
		[&](const char *arg[])
		{
			compare = true;
			return 0;
		}
	);

	po.addOption("-A","--acceptance","acceptance test",1,
		[&](const char *arg[])
		{
			params.acceptTest = true;
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

	if (help || argc<=1)
	{
		std::cerr << "usage: bgmc <options> ..." << std::endl;
		std::cerr << "options:" << std::endl;
		po.printDescriptions(std::cerr);
		return argc<=1 ? 1 : 0;
	}

	cf.addAttribute("N",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGMCParameters*>(pData)->particleCount))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("nmax",
		[](std::istream &is, const char *, void *pData)
		{
			BGMCParameters *pp = static_cast<BGMCParameters*>(pData);
			if (!(is>>pp->nMax))
				return 1;
			pp->beta = kCutoffConst/pp->nMax/pp->nMax;
			return 0;
		}
	);

	cf.addAttribute("beta",
		[](std::istream &is, const char *, void *pData)
		{
			BGMCParameters *pp = static_cast<BGMCParameters*>(pData);
			if (!(is>>pp->beta))
				return 1;
			pp->nMax = static_cast<int32_t>(sqrt(kCutoffConst/pp->beta));
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

	/*std::cerr << "particleCount = " << params.particleCount << std::endl;
	std::cerr << "Nmax = " << params.nMax << std::endl;
	std::cerr << "beta = " << params.beta << std::endl;
	std::cerr << "seed = " << params.seed << std::endl;*/

	if (compare)
	{
		params.sort = false;
		if (r=bgCFCompare(params))
			return r;
	}
	else
	if (r=bgSingleCFSimulation(params))
		return r;

	return 0;
}
