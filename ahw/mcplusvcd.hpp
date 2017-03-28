/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mcplusvcd.hpp
 * Author: posypkin
 *
 * Created on March 11, 2017, 6:39 PM
 */

#ifndef MCPLUSVCD_HPP
#define MCPLUSVCD_HPP

#include <memory>
#include <common/vec.hpp>
#include <funcscale.hpp>
#include <funccnt.hpp>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <methods/varcoorgrad/varcoorgrad.hpp>
#include <pointgen/randpointgen.hpp>
#include <spacefill/spacefillsearch.hpp>

/**
 * Monte-Carlo + Variable Coordinate Descent
 */
class MCplusVCD : public COMPI::Solver<double> {
public:

    class MyWatcher : public BBSEARCH::SpaceFillSearch<double>::Watcher {

        void beforeLocalSearch(double bestf, double inif, int n, const double* x, int cnt) override {
            std::cout << "New point = " << inif << "\n";            
        }

        void update(double prevf, double bestf, int n, const double* prevx, const double* newx, int cnt) override {
            std::cout << "===========================\n";
            std::cout << "Record update from  " << prevf << " to " << bestf << " on step " << cnt << "\n";
            std::cout << "x = " << snowgoose::VecUtils::vecPrint(n, newx) << "\n";
            std::cout << "===========================\n";
        }
    };

    /**
     * Constructor 
     * @param prob problem to process
     * @param minGran minimal granularity for coordinate descent (stopping criterium)
     * @param numPoints number of points for MC search
     */
    MCplusVCD(COMPI::MPProblem<double>& prob, double minGran, int numPoints)
    : mProb(prob), mMinimalGranularity(minGran) {
        const int n = prob.mVarTypes.size();
        std::vector<double> scale(n);
        for (int i = 0; i < n; i++) {
            scale[i] = prob.mBox->mB[i] - prob.mBox->mA[i];
            prob.mBox->mA[i] /= scale[i];
            prob.mBox->mB[i] /= scale[i];
        }
        prob.mObjectives[0] = new COMPI::FunctorScale<double>(*(prob.mObjectives[0]), scale);
        prob.mObjectives[0] = new COMPI::FuncCnt<double>(*(prob.mObjectives[0]));
        std::cout << "new box is " << snowgoose::BoxUtils::toString(*(prob.mBox)) << "\n";

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

#endif
#if 0   
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

        desc->getOptions().mHInit = .5;
        desc->getOptions().mHJ = .2;
#endif
#if 1   
        LOCSEARCH::VarCoorGrad<double>* desc = new LOCSEARCH::VarCoorGrad<double>(prob);


        auto watcher = [&](double fval, const double* x, const std::vector<double>& gran, int stepn) {
            const int n = mProb.mVarTypes.size();
            std::cout << "\n";
            std::cout << "Step: " << stepn << ", ";
            std::cout << "Func calls: " << (dynamic_cast<COMPI::FuncCnt<double>*> (mProb.mObjectives[0]))->mCounters.mFuncCalls << "\n";
            std::cout << "Objective = " << fval << "\n";
            //std::cout << "Solution: " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
            std::cout << "Granularity vector: " << snowgoose::VecUtils::vecPrint(gran.size(), gran.data()) << "\n";
            std::cout << "Maximal granularity: " << snowgoose::VecUtils::maxAbs(n, gran.data(), nullptr) << "\n";
            std::cout.flush();
        };

        desc->getOptions().mHInit = .5;
        desc->getOptions().mHLB = 1e-4;
        //desc->getOptions().mGradStep = .5;
        desc->getOptions().mGradStep = -1;
        desc->getOptions().mGradMaxSteps = 512;
        desc->getOptions().mGradSpeedup = 1;
#endif
        desc->getWatchers().push_back(watcher);
        snowgoose::RandomPointGenerator<double> *rgen = new snowgoose::RandomPointGenerator<double>(*(prob.mBox), numPoints, 1);
        mSFSearch = std::unique_ptr< BBSEARCH::SpaceFillSearch<double> >(new BBSEARCH::SpaceFillSearch<double> (prob, *rgen, *desc));
        mSFSearch->setWatcher(mWatcher);
    }

    bool search(double* x, double& v) override {
        return mSFSearch->search(x, v);
    }

    std::string about() const override {
        return mSFSearch->about();
    }

private:
    const COMPI::MPProblem<double>& mProb;
    double mMinimalGranularity;
    MyWatcher mWatcher;
    std::unique_ptr<BBSEARCH::SpaceFillSearch<double>> mSFSearch;
};

#endif /* MCPLUSVCD_HPP */

