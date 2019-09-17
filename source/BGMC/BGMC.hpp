#ifndef BGMC_HPP_
#define BGMC_HPP_

#include <cinttypes>
#include <string>
#include <random>
#include <iostream>
#include <iomanip>
#include <complex>

static constexpr double kCutoffConst = 0.58;

struct BGMCParameters
{
	uint32_t particleCount;
	int32_t nMax;
	int32_t extraModePairs;
	double beta;
	double gamma;
	uint32_t seed;
	uint32_t batchCount, batchSize;
	uint32_t skip, skipOutput;
	uint32_t firstBatch;
	std::string interactionType;
	std::vector<double> interactionCoefficients;
	std::string output, outputStyle;
	bool useConstDelta;
	bool sort;
	bool acceptTest;
	double delta;
	double betaRatio;
};

struct AlphaRecord
{
	double E;
	double P;
	double Ed, Pd;
	double distance;
	std::vector<std::complex<double> > alpha;
};

struct SimulationInfo
{
	std::vector<double> pAccept;
	double pAcceptMean, pAcceptStdDev, pAcceptMeanStdDev;
	std::vector<std::vector<double> > n;
	std::vector<double> nMean, nStdDev, nMeanStdDev;
	std::vector<double> asymmetry;
	double asymmetryMean, asymmetryStdDev;
	std::vector<double> energy;
	double energyMean, energyStdDev, energyMeanStdDev;
	std::vector<double> momentum;
	double momentumMean, momentumStdDev, momentumMeanStdDev;
	
	inline void append(double accepted, const std::vector<double> &occupation, double A, double E, double P)
	{
		pAccept.push_back(accepted);
		for (uint32_t i=0;i<occupation.size();++i)
			n[i].push_back(occupation[i]);
		asymmetry.push_back(A);
		energy.push_back(E);
		momentum.push_back(P);
	}
	
	template<typename X> void print(std::ostream &os, X batch)
	{
		os << std::setw(8) << batch
			<< " " << std::setw(15) << pAcceptMean
			<< " " << std::setw(15) << nMean[nMean.size()/2] << std::setw(15) << nMeanStdDev[nMeanStdDev.size()/2] << std::setw(15) << nStdDev[nStdDev.size()/2]
			//<< " " << std::setw(15) << asymmetryMean << std::setw(15) << asymmetryStdDev
			<< " " << std::setw(15) << energyMean << std::setw(15) << energyMeanStdDev << std::setw(15) << energyStdDev
			<< " " << std::setw(15) << momentumMean << std::setw(15) << momentumMeanStdDev << std::setw(15) << momentumStdDev
			<< std::endl;
	}
	
	void process();
	void clear(int32_t nMax);
	static void printHead(std::ostream &os);
	void printOccupation(std::ostream &os);
};

class BGMC
{
public:

	struct Energy
	{
		double totalEnergy;
		double kineticEnergy;
		double interactionEnergy;
	};

	virtual double momentum() = 0;
	virtual void energy(Energy &e) = 0;

	virtual uint32_t excitedStatesOccupation() = 0;

	virtual int32_t steps(std::mt19937 &random, uint32_t count) = 0;

	virtual void generate(std::mt19937 &random) = 0;

	virtual void initialize(const BGMCParameters &params) = 0;
	virtual void release() = 0;

	virtual ~BGMC();
};

#endif
