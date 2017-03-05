#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <pointgen/randpointgen.hpp>
#include <memory>
#include "acoustics_bottom_R_weighted.hpp"

class MCplusVCD : public COMPI::Solver<double> {
public:

    MCplusVCD(const COMPI::MPProblem<double>& prob, double minGran, int numPoints)
    : mProb(prob), mMinimalGranularity(minGran) {
        const int n = prob.mVarTypes.size();

#if 0
        LOCSEARCH::CoorDesc<double> desc(*prob, [&](double xdiff, double fdiff, double gran, double fval, int n) {
            return false;
        });

        auto watcher = [&](double xdiff, double fdiff, double gran, double fval, int stp) {
            std::cout << "\n";
            std::cout << "Step: " << stp << ", ";
            std::cout << "Objective = " << fval << "\n";
            std::cout << "Solution: " << << snowgoose::VecUtils::vecPrint(n, x) << "\n";
            std::cout << "Granularity: " << gran << "\n";
            std::cout.flush();
        };


#else    
        LOCSEARCH::VarCoorDesc<double>* desc = new LOCSEARCH::VarCoorDesc<double>(prob, [&](double xdiff, double fdiff, const std::vector<double>& gran, double fval, int n) {
            double a = snowgoose::VecUtils::maxAbs(gran.size(), gran.data());

            std::cout << "maximal granularity = " << a << "\n";
            if (a < mMinimalGranularity)
                return true;
            else
                return false;
            });

        auto watcher = [&](double xdiff, double fdiff, const std::vector<double>& gran, double fval, int stp) {
            const int n = mProb.mVarTypes.size();
            std::cout << "\n";
            std::cout << "Step: " << stp << ", ";
            std::cout << "Objective = " << fval << "\n";
            //std::cout << "Solution: " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
            std::cout << "Granularity vector: " << snowgoose::VecUtils::vecPrint(gran.size(), gran.data()) << "\n";
            std::cout.flush();
        };

#endif
        desc->getWatchers().push_back(watcher);

        for (int i = 0; i < n; i++) {
            desc->getOptions().mShifts[i] = prob.mBox->mB[i] - prob.mBox->mA[i];
        }
        mLocalSearch = std::move(std::unique_ptr<COMPI::Solver<double>>(desc));
        snowgoose::RandomPointGenerator<double> *rgen = new snowgoose::RandomPointGenerator<double>(*(prob.mBox), numPoints, 1);
        mPointGenerator = std::move(std::unique_ptr<snowgoose::PointGenerator<double>>(rgen));
    }

    bool search(double* x, double& v) override {
        const int n = mProb.mVarTypes.size();
        double* tx = new double[n];
        while (mPointGenerator->getPoint(tx)) {

            double tv = mProb.mObjectives.at(0)->func(tx);
            std::cout << "Initial value: " << tv << "\n";
            std::cout << "Initial Solution: " << snowgoose::VecUtils::vecPrint(n, tx) << "\n";
            std::cout.flush();
            mLocalSearch->search(tx, tv);
            std::cout << "Found value: " << tv << "\n";
            std::cout << "Solution found: " << snowgoose::VecUtils::vecPrint(n, tx) << "\n";
            if (tv < v) {
                snowgoose::VecUtils::vecCopy(n, tx, x);
                v = tv;
            }
        }
        delete [] tx;
    }

    std::string about() const override {
        return "Monte-Carlo + " + mLocalSearch->about();
    }

private:
    const COMPI::MPProblem<double>& mProb;
    std::unique_ptr<COMPI::Solver<double>> mLocalSearch;
    std::unique_ptr<snowgoose::PointGenerator<double>> mPointGenerator;
    double mMinimalGranularity;
};

int main(int argc, char *argv[]) {
    // Setup problem
    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.verbosity = 0;
    sspemdd_seq.readScenario("312_bottom_R_weighted260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterProblemFactory ahwpf(vPair);
    for (auto p : vPair) {
        std::cout << p.first << " : " << p.second << "\n";
    }
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();

    // Setup solver
    const int n = 3;
    const int numberOfPoints = 4;
    const double minimalGranularity = 1e-4;
    double x[n];
    snowgoose::BoxUtils::getCenter(*(prob->mBox), (double*) x);
    double v = prob->mObjectives.at(0)->func(x);

    MCplusVCD mcsearch(*prob, minimalGranularity, numberOfPoints);

    // Run solver
    std::cout << "Searching with " << mcsearch.about() << "\n";
    mcsearch.search(x, v);
    
    // Print results
    std::cout << "Found x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Objective value = " << v << "\n";
    return 0;
}