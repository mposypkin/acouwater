#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <pointgen/randpointgen.hpp>
#include <spacefill/spacefillsearch.hpp>
#include "mcplusvcd.hpp"
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
    COMPI::MPProblem<double>* prob = ahwupf.getProblem();

    // Setup solver
    const int n = prob->mVarTypes.size();
    const int numberOfPoints = 4;
    const double minimalGranularity = 1e-4;

    MCplusVCD mcsearch(*prob, minimalGranularity, numberOfPoints);

    // Setup initial point
    double x[n];
    snowgoose::BoxUtils::getCenter(*(prob->mBox), (double*) x);
    double v = prob->mObjectives.at(0)->func(x);

    // Run solver
    std::cout << "Searching with " << mcsearch.about() << "\n";
    mcsearch.search(x, v);

    // Print results
    std::cout << "Found x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Objective value = " << v << "\n";
    return 0;


    return 0;
}