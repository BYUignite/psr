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
    psr.setT(2220.42);

    //--------------- solve

    psr.setInitialGuess();
    psr.solve();

    cout << endl << psr.tau << endl;
    cout << endl << psr.y[gas->speciesIndex("CO2")] << endl;
    cout << endl << psr.y[gas->speciesIndex("CO")] << endl;

    return 0;
}
