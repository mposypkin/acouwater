/* 
 * File:   acoustics_homog_water.hpp
 * Author: Oleg Zaikin
 *
 * Created on Devember 9, 2016, 8:12 PM
 */

#ifndef ACOUSTICS_HOMOG_WATER_HPP
#define ACOUSTICS_HOMOG_WATER_HPP

#include <mpproblem.hpp>
#include <box/box.hpp>
#include <vector>
#include <sspemdd_sequential.h>

namespace ACOUSTIC {

    class AcousticsHomogWaterObjective : public COMPI::Functor <double> {
    public:

        AcousticsHomogWaterObjective(int n) : mN(n) {
        }

        double func(const double* x) {
			sspemdd_sequential sspemdd_seq;
			sspemdd_seq.readScenario("34_bottom_R_weighted160.txt");
			sspemdd_seq.readInputDataFromFiles();
            sspemdd_seq.init();
            search_space_point cur_point;
			// set const values
			cur_point.tau = sspemdd_seq.tau1;
			cur_point.cws = sspemdd_seq.cw1_arr;
			// variable values
			cur_point.R = x[0];
            cur_point.rhob = x[1];
            cur_point.cb = x[2];
            return sspemdd_seq.fill_data_compute_residual(cur_point);
        }

    private:
        int mN;
    };

    /**
     * Factory to produce instances of DeJong optimization problem
     */
    class AcousticsHomogWaterProblemFactory {
    public:

        /**
         * Constructor
         * @param n problem dimension
         * @param vecA - a vector of left borders for each dimension
         * @param vecB - a vector of right borders for each dimension
         * n = 3
         * x0 = cb in (1600, 1900)
         * x1 = rhob in (1.4, 2.0)
         * x2 = tau in (-0.015, 0.015)
         * vecA = {1600, 1.4, -0.015}
         * vecB = {1900, 2.0, 0.015}
         */
        AcousticsHomogWaterProblemFactory(int n, const std::vector<double>& vecA, const std::vector<double>& vecB) :
        mN(n),
        mA(vecA),
        mB(vecB) {
        }

        COMPI::MPProblem<double>* getProblem() const {
            COMPI::MPProblem<double>* prob = new COMPI::MPProblem<double>();
            prob->mVarTypes.assign(mN, COMPI::MPProblem<double>::VariableTypes::GENERIC);
            prob->mObjectives.push_back(new AcousticsHomogWaterObjective(mN));
            prob->mBox = new snowgoose::Box<double>(mN);
            for (int i = 0; i < mN; i++) {
                prob->mBox->mA[i] = mA[i];
                prob->mBox->mB[i] = mB[i];
            }
            return prob;
        }

    private:
        int mN;
        std::vector<double> mA;
        std::vector<double> mB;
    };
}

#endif /* ACOUSTICS_HOMOG_WATER_HPP */

