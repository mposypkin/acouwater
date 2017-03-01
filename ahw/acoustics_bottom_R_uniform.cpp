#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <pointgen/randpointgen.hpp>
#include "acoustics_bottom_R_uniform.hpp"

int main(int argc, char *argv[]) {
    const int n = 3;
    const int numberOfPoints = 4;
    const double minimalGranularity = 1e-4;

    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.verbosity = 0;
    sspemdd_seq.readScenario("311_bottom_R_uniform260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterUniformProblemFactory ahwupf(vPair);
    COMPI::MPProblem<double>* prob = ahwupf.getProblem();

    const double eps = 1e-5;
    const double vref = 0.0130504;
    double x[n] = {7016, 1.7, 1735};
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";
    if (SGABS(v - vref) >= eps)
        return -1;

    LOCSEARCH::VarCoorDesc<double> desc(*prob, [&](double xdiff, double fdiff, const std::vector<double>& gran, double fval, int n) {
        double a = snowgoose::VecUtils::maxAbs(gran.size(), gran.data());

        std::cout << "maximal granularity = " << a << "\n";
        if (a < minimalGranularity)
            return true;
        else
            return false;
        });

    auto watcher = [&](double xdiff, double fdiff, const std::vector<double>& gran, double fval, int stp) {
        std::cout << "\n";
        std::cout << "Step: " << stp << ", ";
        std::cout << "Objective = " << fval << "\n";
        std::cout << "Solution: " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        std::cout << "Granularity vector: " << snowgoose::VecUtils::vecPrint(gran.size(), gran.data()) << "\n";
        std::cout.flush();
    };

    desc.getWatchers().push_back(watcher);


    for (int i = 0; i < n; i++) {
        desc.getOptions().mShifts[i] = vPair[i].second - vPair[i].first;
    }
    
    snowgoose::RandomPointGenerator<double> rgen(*(prob->mBox), 1);
    for (int i = 0; i < numberOfPoints; i++) {
        rgen.getPoint(x);

        double v = prob->mObjectives.at(0)->func(x);
        std::cout << "Initial value: " << v << "\n";
        std::cout << "Initial Solution: " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        std::cout.flush();

        desc.search(x, v);
        std::cout << "Found value: " << v << "\n";
        std::cout << "Solution found: " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    }


    return 0;
}