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
	params.batchCount = 9;
	params.batchSize = 1024*64;
	params.skip = 4;
	params.interactionType = "contact";

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

	po.addOption("-n","--nmax <count>","set Nmax",2,
		[&](const char *arg[])
		{
			params.nMax = std::stol(arg[1]);
			params.beta = kCutoffConst/params.nMax/params.nMax;
			return 0;
		}
	);

	po.addOption("-e","--extra <count>","extra mode pairs",2,
		[&](const char *arg[])
		{
			params.extraModePairs = std::stol(arg[1]);
			return 0;
		}
	);

	po.addOption("-b","--beta <number>","set beta(inverse temperature)",2,
		[&](const char *arg[])
		{
			params.beta = std::stod(arg[1]);
			params.nMax = static_cast<int32_t>(sqrt(kCutoffConst/params.beta));
			return 0;
		}
	);

	po.addOption("-T","--temperature <number>","set temperature)",2,
		[&](const char *arg[])
		{
			params.beta = 1.0/std::stod(arg[1]);
			params.nMax = static_cast<int32_t>(sqrt(kCutoffConst/params.beta));
			return 0;
		}
	);

	po.addOption("-g","--gamma <number>","interaction strength",2,
		[&](const char *arg[])
		{
			params.gamma = std::stod(arg[1]);
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
		std::cerr << "usage: bgmc <options> ..." << std::endl;
		std::cerr << "options:" << std::endl;
		po.printDescriptions(std::cerr);
		return 0;
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

	if (r=bgSimulationCF(params))
		return r;

	return 0;
}
