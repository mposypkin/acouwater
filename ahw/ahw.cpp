#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include "acoustics_homog_water.hpp"

int main(int argc, char *argv[]) {
    const int n = 3;
    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.readScenario("312_bottom_R_weighted260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterProblemFactory ahwpf(vPair);
    for(auto p : vPair) {
        std::cout << p.first << " : " <<  p.second << "\n";
    }
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();

    const double eps = 1e-10;
    const double vref = 3.85835e-07;
    //double x[n] = {7019, 1.75, 1705};
    double x[n];
    x[0]= 0.5 * (sspemdd_seq.R1 + sspemdd_seq.R2);
    x[1] = 0.5 * (sspemdd_seq.rhob1 + sspemdd_seq.rhob2);
    x[2] = 0.5 * (sspemdd_seq.cb1 + sspemdd_seq.cb2);
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";
    std::cout.flush();

  
    int cnt = 0;
    
   
#if 0
    LOCSEARCH::CoorDesc<double> desc(*prob, [&](double xdiff, double fdiff, double gran, double fval, int n) {
        cnt++; 
        std::cout << "CNT = " << cnt << ", ";
        std::cout << "VAL = " << fval << "\n";
        std::cout.flush();
        return false;
    });
#else    
     LOCSEARCH::VarCoorDesc<double> desc(*prob, [&](double xdiff, double fdiff, const std::vector<double>& gran, double fval, int n) {
        cnt++; 
        std::cout << "CNT = " << cnt << ", ";
        std::cout << "VAL = " << fval << "\n";
        std::cout.flush();
        return false;
    });
#endif
    
    for(int i = 0; i < n; i ++) {
        desc.getOptions().mShifts[i] = vPair[i].second - vPair[i].first;
    }

    desc.search(x, v);
    std::cout << "v = " << v << "\n";
    
    return 0;
}