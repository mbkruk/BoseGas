#include "BGMC.hpp"
#include "../BGCommon/Math.hpp"

void SimulationInfo::process()
{
	meanStdDev(pAccept,&pAcceptMean,&pAcceptStdDev);
	pAcceptMeanStdDev = pAcceptStdDev/sqrt(pAccept.size()-1);
	for (uint_fast32_t i=0;i<n.size();++i)
	{
		meanStdDev(n[i],&nMean[i],&nStdDev[i]);
		nMeanStdDev[i] = nStdDev[i]/sqrt(n[i].size()-1);
	}
	meanStdDev(asymmetry,&asymmetryMean,&asymmetryStdDev);
	meanStdDev(energy,&energyMean,&energyStdDev);
	energyMeanStdDev = energyStdDev/sqrt(energy.size()-1);
	meanStdDev(momentum,&momentumMean,&momentumStdDev);
	momentumMeanStdDev = momentumStdDev/sqrt(momentum.size()-1);
}

void SimulationInfo::clear(int32_t nMax)
{
	pAccept.clear();
	n.resize(2*nMax+1);
	nMean.resize(n.size());
	nStdDev.resize(n.size());
	nMeanStdDev.resize(n.size());
	for (uint_fast32_t i=0;i<n.size();++i)
	if (!n[i].empty())
		n[i].clear();
	asymmetry.clear();
	energy.clear();
	momentum.clear();
}

void SimulationInfo::printHead(std::ostream &os)
{
	os << std::setw(8) << "batch"
		<< " " << std::setw(15) << "P accept"
		<< " " << std::setw(15) << "n0 mean" << std::setw(15) << "+/-" << std::setw(15) << "n0 std dev"
		//<< " " << std::setw(15) << "asym mean" << std::setw(15) << "asym std dev"
		<< " " << std::setw(15) << "E mean" << std::setw(15) << "+/-" << std::setw(15) << "E std dev"
		<< " " << std::setw(15) << "P mean" << std::setw(15) << "+/-" << std::setw(15) << "P std dev"
		<< std::endl;
}

void SimulationInfo::printOccupation(std::ostream &os)
{
	os << std::setw(8) << "mode"
		<< " " << std::setw(15) << "occupation"
		<< " " << std::setw(15) << "+/-"
		<< " " << std::setw(15) << "std dev"
		<< std::endl;

	for (uint_fast32_t i=0;i<n.size();++i)
	{
		os << std::setw(8) << ((int32_t)i-(int32_t)n.size()/2)
			<< " " << std::setw(15) << nMean[i]
			<< " " << std::setw(15) << nMeanStdDev[i]
			<< " " << std::setw(15) << nStdDev[i]
			<< std::endl;
	}
}

BGMC::~BGMC()
{
}
