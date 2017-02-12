#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include "acoustics_homog_water.hpp"

int main(int argc, char *argv[])
{
	const int n = 3;
	sspemdd_sequential sspemdd_seq;
	sspemdd_seq.readScenario("312_bottom_R_weighted260.txt");
	std::vector<std::pair<double, double>> vPair;
	vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
	vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
	vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
	ACOUSTIC::AcousticsHomogWaterProblemFactory ahwpf(vPair);
	COMPI::MPProblem<double>* prob = ahwpf.getProblem();

	const double eps = 1e-10;
	const double vref = 3.85835e-07;
	double x[n] = { 7019, 1.75, 1705 };
	double v = prob->mObjectives.at(0)->func(x);
	std::cout << "v = " << v << "\n";
	if ( SGABS(v - vref) >= eps)
		return -1;

	return 0;
}