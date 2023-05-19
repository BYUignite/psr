#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"

#include "kinsol_wrapper.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class PSR {

    //////////////////// DATA MEMBERS //////////////////////

public:

    std::shared_ptr<Cantera::ThermoPhase> gas;
    std::shared_ptr<Cantera::Kinetics>    kin;

    std::vector<double> yin;
    double              hin;
    double              P;
    double              T;

    size_t neq;

    std::vector<double> y;
    double tau;

    kinsol_wrapper *kwrap;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    PSR(std::shared_ptr<Cantera::ThermoPhase> gas_, 
            std::shared_ptr<Cantera::Kinetics> kin_);

    //////////////////// MEMBER FUNCTIONS /////////////////

    void setInlet(std::vector<double> &yin_, double hin_, double P_) {
        yin = yin_; 
        hin = hin_;
        P = P_; 
    }

    void setT(double T_) { T = T_; }
    void setInitialGuess();
    void setInitialGuess(const double taug, const std::vector<double> yg, const bool doEq=false);
    void solve();
    void Func(const double *ytau, double *F);
};
