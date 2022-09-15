#include "psr.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

int main() {

    //--------------- user inputs

    string mechName   = "gri30.yaml";

    string xin = "CH4:1, O2:2, N2:7.52";
    double Tin = 300;
    double P   = 101325;

    int    nT   = 500;            // number of T values to solve for

    double Tmin = Tin;
    double TmaxDelta = -5.1;      // this can be 0.1 or 0.01 for stoich methane/air, but higher like 5 or more for lean to 0.03 mixf

    double taug = 0.01;           // first guess for tau; solver likes a low guess
        
    //--------------- initialize cantera

    auto sol = newSolution(mechName, "", "None");
    auto gas = sol->thermo();
    auto kin = sol->kinetics();

    size_t nsp  = gas->nSpecies();

    //--------------- inlet gas state

    vector<double> yin(nsp);
    double         hin;

    gas->setState_TPX(Tin, P, xin);
    gas->getMassFractions(&yin[0]);
    hin = gas->enthalpy_mass();

    //--------------- adiabatic equilibrium gas state

    gas->equilibrate("HP");
    double Tad = gas->temperature();
    vector<double> yad(nsp);
    gas->getMassFractions(&yad[0]);

    //--------------- PSR object, scaling arrays

    PSR psr(gas, kin, yin, hin, P);

    vector<double> y_tau_scale(nsp+1, 1.0);
    vector<double> f_scale(nsp+1, 1.0);

    //for(int k=0; k<nvar; k++){{{{
    //    ytauscale[k] = ytau[k];
    //    if(ytauscale[k] < 1E-10) ytauscale[k] = 1E-10;
    //    ytauscale[k] = 1.0/ytauscale[k];
    //    //fscale[k] = 1.0/(0.1/tau);
    //}}}}

    //--------------- solve psr

    double Tmax = Tad + TmaxDelta;
    vector<double> Tvec(nT);      // temperature values
    for(int i=0; i<nT; i++)
        Tvec[i] = Tmax - (double)(i)/(nT-1) * (Tmax - Tmin);

    vector<double> y_tau = yad;   // unknown vector: species mass fractions and tau
    y_tau.push_back(taug);

    vector<vector<double>> ys(nT, vector<double>(nsp, 0.0));   // store y vec at each point
    vector<double>         taus(nT);                           // store tau at each point

    ofstream ofile("tau_T.dat");
    ofile << "# tau (s), T (K) " << endl;

    for(int i=0; i<nT; i++) {                                      // loop over each point
        psr.setT(Tvec[i]);
        psr.solvePSR(y_tau, y_tau_scale, f_scale);                 // solve psr at this point
        ys[i].insert(ys[i].begin(), y_tau.begin(), y_tau.end()-1); // store solution
        taus[i] = y_tau.back();                                    // store solution
        ofile << y_tau.back() << " " << Tvec[i] << endl;
    }

    ofile.close();
    system("ipython3 plot.py");

    return 0;
}

