#include <iostream>
#include <sspemdd_sequential.h>
#include <pointgen/randpointgen.hpp>
#include <methods/mtcoordescent/ctcoordescent.hpp>
#include <box/boxutils.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "acoustics_bottom_R_uniform.hpp"

int main(int argc, char *argv[]) {
    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.verbosity = 0;
    sspemdd_seq.readScenario("311_bottom_R_uniform260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterUniformProblemFactory ahwupf(vPair);
	
	COMPI::MPProblem<double> *mpp = ahwupf.getProblem();

	LOCSEARCH::CTCoordinateDescent<double> desc(*mpp);
	desc.getOptions().mHInit = .1;
	//desc.getOptions().mHLB = 1e-10;
	//desc.getOptions().mParallelMode = false;
	desc.getOptions().mNumThreads = 2;
	desc.getOptions().mGradLB = 0;
    
    // Setup initial point
	const int n = mpp->mVarTypes.size();
	std::cout << "n " << n << "\n";
	double x[n];
	snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);

	// Search
	double v;
	auto start = std::chrono::steady_clock::now();
	std::cout << "x : ";
	for (int i = 0; i < n; i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
	std::cout << "start search" << std::endl;
	bool rv = desc.search(x, v);
	auto diff = std::chrono::steady_clock::now() - start;
	std::cout << desc.about() << "\n";
	std::cout << "Found v = " << mpp->mObjectives.at(0)->func(x) << "\n";
	std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
	//std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
	std::cout << "Execution time = " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;

    // Rescale x
    /*const double *s = mcsearch.getScale().data();
    snowgoose::VecUtils::vecMultVect(n, x, s, x);
    // Print results
    std::cout << "Found x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Objective value = " << v << "\n";*/

    return 0;
}