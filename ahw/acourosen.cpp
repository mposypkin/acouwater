#include <iostream>
#include <utility>
#include <chrono>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <methods/rosenbrock/rosenbrockmethod.hpp>
#include "acoustics_bottom_R_uniform.hpp"

double inc;
double dec;
bool doort;

void rosenSearch(double& v, double* x, const COMPI::MPProblem<double>& prob) {
    LOCSEARCH::RosenbrockMethod<double> desc(prob);
    desc.getOptions().mHInit = std::vector<double>(prob.mBox->mDim, 1.);
    desc.getOptions().mMaxStepsNumber = 10000;
    desc.getOptions().mMinGrad = 1e-3;
    desc.getOptions().mHLB = 1e-4 * desc.getOptions().mMinGrad;
    desc.getOptions().mDec = dec;
    desc.getOptions().mInc = inc;
    desc.getOptions().mDoOrt = doort;
    //    desc.getOptions().mDoOrt = false;
    desc.getOptions().mDoTracing = true;
    desc.getWatchers().push_back([](double fval, const double* x, const std::vector<double>& gran, bool success, double grad, double* dirs, int stpn) {
        std::cout << "stage " << stpn << ", fval = " << fval << std::endl;
        std::cout << "gran = (";
        for (auto g : gran) {
            std::cout << " " << g;
        }
        std::cout << " )" << std::endl;
        if (success) {
            std::cout << "successful stage, grad = " << grad << std::endl;
        } else {
            std::cout << "unsuccessful stage, dirs =" << std::endl;
                    const int n = gran.size();
            for (int i = 0; i < n; i++) {
                const double * const p = dirs + i * n;
                for (int j = 0; j < n; j++) {

                    std::cout << " " << p[j];
                }
                std::cout << std::endl;
            }
        }
        std::cout << "-----" << std::endl;
    });
    std::cout << desc.about();
    std::cout << "before v = " << v << std::endl;
    bool rv = desc.search(x, v);
    std::cout << "after v = " << v << std::endl;
}

int main(int argc, char *argv[]) {
    constexpr int n = 3;

    dec = atof(argv[1]);
    inc = atof(argv[2]);
    doort = (atoi(argv[3]) == 0) ? false : true;

    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.verbosity = 0;
    sspemdd_seq.readScenario("311_bottom_R_uniform260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterUniformProblemFactory ahwpf(vPair);
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();


    double x[n];
    snowgoose::BoxUtils::getCenter(*(prob->mBox), x);
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";

    auto start = std::chrono::steady_clock::now();
    rosenSearch(v, x, *prob);
    auto diff = std::chrono::steady_clock::now() - start;
    std::cout << "Execution time = " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;

    return 0;
}