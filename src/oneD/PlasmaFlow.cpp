//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/PlasmaFlow.h"

using namespace std;

namespace Cantera
{

PlasmaFlow::PlasmaFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    IonFlow(ph, nsp, points)
{
    vector<string> list;
    list.push_back("N2");
    list.push_back("N");
    list.push_back("O2");
    list.push_back("O");
    list.push_back("H2");
    list.push_back("H");
    list.push_back("CO2");
    list.push_back("CO");
    list.push_back("H2O");
    list.push_back("CH4");
    // check valid index 
    for (size_t i = 0; i < list.size(); i++) {
        size_t k = m_thermo->speciesIndex(list[i]);
        if (k != npos) {
            m_collisionSpeciesIndex.push_back(k);
        }
    }
    zdplaskin_init();
    //zdplaskin_set_density("O2", 0.0000001);
    //cout << getElectronTemperature() << endl;
}

void PlasmaFlow::resize(size_t components, size_t points){
    IonFlow::resize(components, points);
}

void PlasmaFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x, j0, j1);
}

void PlasmaFlow::eval(size_t jg, double* xg,
                  double* rg, integer* diagg, double rdt)
{
    IonFlow::eval(jg, xg, rg, diagg, rdt);
     // start of local part of global arrays
    double* x = xg + loc();
    //double* rsd = rg + loc();
    //integer* diag = diagg + loc();
    // Define boundary Indexes
    size_t jmin, jmax, j0, j1;
    getBoundaryIndexes(jg, jmin, jmax, j0, j1);
    for (size_t j = jmin; j <= jmax; j++) {
        zdplaskin_set_conditions(T(x, j), 0.0);
        for (size_t k : m_collisionSpeciesIndex) {
            zdplaskin_set_density(m_thermo->speciesName(k), X(x, k, j));
        }
    }
}

void PlasmaFlow::evalEEDF(size_t j, double* x, double* rsd, integer* diag, double rdt)
{       

}

}
