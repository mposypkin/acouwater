#include <iostream>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include "acoustics_homog_water.hpp"

int main(int argc, char *argv[])
{
    const int n = 3;
    const double eps = 0.001;
    const double vref = 100;
    std::vector<double> a;
    std::vector<double> b;
    a.assign(n, 1);
    b.assign(n, 2);
    ACOUSTIC::AcousticsHomogWaterProblemFactory ahwpf(n, a, b);
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();
    double x[n] = {1, 1, 1};
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";
    SG_ASSERT(SGABS(v - vref) < eps);
    return 0;
}