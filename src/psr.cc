#include "psr.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////

PSR::PSR(std::shared_ptr<Cantera::ThermoPhase> gas_, 
         std::shared_ptr<Cantera::Kinetics>    kin_) :
    gas(gas_),
    kin(kin_) {

    neq = gas->nSpecies()+1;

    yin.resize(gas->nSpecies());
    y.resize(gas->nSpecies());

    //----------

    vector<double> scales_v(neq, 1.0);
    vector<double> scales_f(neq, 1.0);
    vector<double> constraints(neq, 0.0); constraints.back() = 1.0;

    kwrap = new kinsol_wrapper(this, neq, scales_v, scales_f, constraints);
}

////////////////////////////////////////////////////////////////////////////////

void PSR::setInitialGuess(const double taug, const vector<double> yg, const bool doEq) {

    tau = taug;
    if (doEq) {
        gas->setMassFractions(&yin[0]);
        gas->setState_HP(hin, P);
        gas->equilibrate("HP");
        gas->getMassFractions(&y[0]);
    }
    else {
        y = yg;
    }
}

////////////////////////////////////////////////////////////////////////////////

void PSR::solve() {
    
    vector<double> ytau = y;
    ytau.push_back(tau);

    kwrap->solve(ytau);

    for(int i=0; i<neq-1; i++)
        y[i] = ytau[i];
    tau = ytau.back();
}

////////////////////////////////////////////////////////////////////////////////

void PSR::Func(const double *ytau, double *F) {

    gas->setMassFractions_NoNorm(ytau);
    gas->setState_TP(T, P);
    double rho = gas->density();
    vector<double> rr(neq-1);
    kin->getNetProductionRates(&rr[0]);

    for(size_t k=0; k<neq-1; k++)
        F[k] = (yin[k] - ytau[k])/ytau[neq-1] + rr[k]*gas->molecularWeight(k)/rho;
    F[neq-1] = gas->enthalpy_mass() - hin;

}

////////////////////////////////////////////////////////////////////////////////
// Kinsol interface; Kinsol calls this function, which then calls user_data's Func 

int FuncKinsol(N_Vector yvec, N_Vector fvec, void *user_data) {

    PSR *psr = static_cast<PSR *>(user_data);

    double *y = N_VGetArrayPointer(yvec);
    double *F = N_VGetArrayPointer(fvec);

    psr->Func(y, F);

    return 0;
}
