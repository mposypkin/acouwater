#include <iostream>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include "acoustics_homog_water.hpp"

int main(int argc, char *argv[])
{
    const int n = 3;
    const double eps = 1e-10;
    const double vref = 3.85835e-07;
    std::vector<double> a;
    std::vector<double> b;
	sspemdd_sequential sspemdd_seq;
	sspemdd_seq.readScenario("312_bottom_R_weighted260.txt");
	a.resize(n);
	a[0] = sspemdd_seq.R1;
	a[1] = sspemdd_seq.rhob1;
	a[2] = sspemdd_seq.cb1;
	b.resize(n);
	b[0] = sspemdd_seq.R2;
	b[1] = sspemdd_seq.rhob2;
	b[2] = sspemdd_seq.cb2;
    ACOUSTIC::AcousticsHomogWaterProblemFactory ahwpf(n, a, b);
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();
    double x[n] = {7019, 1.75, 1705};
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";
    SG_ASSERT(SGABS(v - vref) < eps);
    return 0;
}