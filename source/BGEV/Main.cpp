#include <iostream>
#include <cinttypes>

#include "../BGCommon/ProgramOptions.hpp"
#include "../BGCommon/ConfigFile.hpp"
#include "BGEV.hpp"

int main(int argc, const char *argv[])
{
	ProgramOptions po;
	ConfigFile cf;
	int r;
	bool help = false;
	std::vector<std::string> inputFiles;
	BGEVParameters params;

	params.particleCount = 1;
	params.nMax = 0;
	params.gamma = 1.0;

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
			return 0;
		}
	);

	po.addOption("-g","--gamma <number>","set gamma(interaction strength)",2,
		[&](const char *arg[])
		{
			params.gamma = std::stod(arg[1]);
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

	cf.addAttribute("N",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->particleCount))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("nmax",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->nMax))
				return 1;
			return 0;
		}
	);

	cf.addAttribute("gamma",
		[](std::istream &is, const char *, void *pData)
		{
			if (!(is>>static_cast<BGEVParameters*>(pData)->gamma))
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

	std::cerr << "particleCount = " << params.particleCount << std::endl;
	std::cerr << "Nmax = " << params.nMax << std::endl;
	std::cerr << "gamma = " << params.gamma << std::endl;

	return 0;
}
