#include "psr.h"

#include <string>
#include <iostream>
#include <iomanip>
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
    double TmaxDelta = -5.1;      // this can be 0.1 or 0.01 for stoich methane/air, but higher like 5 or more for lean to 0.03 mixf

    double taug = 0.01;           // first guess for tau; solver likes a low guess

    //--------------- initialize cantera

    auto sol = newSolution(mechName, "", "none");
    auto gas = sol->thermo();
    auto kin = sol->kinetics();

    size_t nsp  = gas->nSpecies();

    //--------------- inlet gas state

    vector<double> yin(nsp);
    double         hin;

    gas->setState_TPX(Tin, P, xin);
    gas->getMassFractions(&yin[0]);
    hin = gas->enthalpy_mass();

    //--------------- psr object

    PSR psr(gas, kin);

    psr.setInlet(yin, hin, P);

    //--------------- adiabatic equilibrium state

    gas->equilibrate("HP");
    double Teq = gas->temperature();

    //--------------- compute S-curve

    double Tmax = Teq + TmaxDelta;
    vector<double> Tvec(nT);      // temperature values
    double Tmin = Tin;
    for(int i=0; i<nT; i++)
        Tvec[i] = Tmax - (double)(i)/(nT-1) * (Tmax - Tmin);

    vector<vector<double>> ys(nT, vector<double>(nsp, 0.0));   // store y vec at each point
    vector<double>         taus(nT);                           // store tau at each point

    ofstream ofile("tau_T.dat");
    ofile << "# tau (s), T (K), Y_CO, Y_CO2 " << endl;

    ofile << scientific;
    ofile << setprecision(10);
    for(int i=0; i<nT; i++) {        // loop over each point
        psr.setT(Tvec[i]);
        i==0 ? psr.setInitialGuess(taug, yin, true) :        // true --> equilibrium not yin
               psr.setInitialGuess(psr.tau, psr.y, false);
        psr.solve();                 // solve psr at this point
        ys[i] = psr.y;               // store solution
        taus[i] = psr.tau;
        ofile << taus[i] << " " 
              << Tvec[i] << " " 
              << ys[i][gas->speciesIndex("CO")] << " "
              << ys[i][gas->speciesIndex("CO2")] << " "
              << endl;
        cout << taus[i] << " " 
              << Tvec[i] << " " 
              << ys[i][gas->speciesIndex("CO")] << " "
              << ys[i][gas->speciesIndex("CO2")] << " "
              << endl;
    }

    ofile.close();
    system("ipython plot.py");

    return 0;
}
