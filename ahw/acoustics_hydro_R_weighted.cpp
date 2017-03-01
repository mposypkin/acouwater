#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <pointgen/randpointgen.hpp>
#include "acoustics_hydro_R_weighted.hpp"

int main(int argc, char *argv[]) {
	const int n = 5;
	sspemdd_sequential sspemdd_seq;
	sspemdd_seq.verbosity = 0;
	sspemdd_seq.readScenario("39_hydro_R_weighted260.txt");
	std::vector<std::pair<double, double>> vPair;
	vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
	vPair.push_back(std::make_pair(sspemdd_seq.cw1_arr[1], sspemdd_seq.cw2_arr[1]));
	vPair.push_back(std::make_pair(sspemdd_seq.cw1_arr[2], sspemdd_seq.cw2_arr[2]));
	vPair.push_back(std::make_pair(sspemdd_seq.cw1_arr[3], sspemdd_seq.cw2_arr[3]));
	vPair.push_back(std::make_pair(sspemdd_seq.cw1_arr[4], sspemdd_seq.cw2_arr[4]));
	
	ACOUSTIC::AcousticsHydroRWeightedProblemFactory ahrwpf(vPair);
	COMPI::MPProblem<double>* prob = ahrwpf.getProblem();

	const double eps = 1e-10;
	const double vref = 3.15003e-07;
	double x[n] = { 6981, 1480, 1485, 1450, 1470 };
	double v = prob->mObjectives.at(0)->func(x);
	std::cout << "v = " << v << "\n";
	if (SGABS(v - vref) >= eps)
		return -1;
}