//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/PlasmaFlow.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/zdplaskin.h"

using namespace std;

namespace Cantera
{

PlasmaFlow::PlasmaFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    IonFlow(ph, nsp, points),
    m_do_plasma(false)
{
    vector<string> list;
    list.push_back("N2");
    list.push_back("O2");
    list.push_back("H2");
    list.push_back("H2O");
    list.push_back("CO2");
    list.push_back("CO");
    list.push_back("CH4");
    // check valid index 
    for (size_t i = 0; i < list.size(); i++) {
        size_t k = m_thermo->speciesIndex(list[i]);
        if (k != npos) {
            m_collisionSpeciesIndex.push_back(k);
        }
    }

    zdplaskinInit();

    for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
        char* cstring[20];
        zdplaskinGetSpeciesName(cstring, &i);
        string speciesName(*cstring);
        size_t k = m_thermo->speciesIndex(speciesName);
        m_plasmaSpeciesIndex.push_back(k);
    }
}

void PlasmaFlow::resize(size_t components, size_t points) {
    IonFlow::resize(components, points);
}

void PlasmaFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x, j0, j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_overwrite_eTransport && (m_kElectron != npos)) {
            if (m_import_electron_transport) {
                updateEEDF(x,j);
                const double number_density = ND_t(x,j);
                m_mobility[m_kElectron+m_nsp*j] = zdplaskinGetElecMobility(&number_density);
                m_diff[m_kElectron+m_nsp*j] = zdplaskinGetElecDiffCoeff();
            } else {
                m_mobility[m_kElectron+m_nsp*j] = 0.4;
                m_diff[m_kElectron+m_nsp*j] = 0.4*(Boltzmann * T(x,j)) / ElectronCharge;
            }
        }
    }
}

void PlasmaFlow::eval(size_t jg, double* xg,
                  double* rg, integer* diagg, double rdt)
{
    IonFlow::eval(jg, xg, rg, diagg, rdt);
    /*
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
    */
}

void PlasmaFlow::solvePlasma()
{
    bool changed = false;
    if (m_do_plasma) {
        changed = true;
    }
    m_do_plasma = true;
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
    }
}

void PlasmaFlow::getWdot(doublereal* x, size_t j)
{
    setGas(x,j);
    m_kin->getNetProductionRates(&m_wdot(0,j));
    if (m_do_plasma) {
        updateEEDF(x, j);
        zdplaskinGetPlasmaSource(&m_wdot_plasma);
        size_t i = 0;
        for (size_t k : m_plasmaSpeciesIndex) {
            i++;
            m_wdot(k,j) += m_wdot_plasma[i];
        }  
    }
}

void PlasmaFlow::updateEEDF(double* x, size_t j)
{
    for (size_t k : m_collisionSpeciesIndex) {
        const char* species[10] = {m_thermo->speciesName(k).c_str()};
        const double num_density = ND(x,k,j);
        zdplaskinSetDensity(species, &num_density);
    }
    const double temperature = T(x,j);
    const double ruduced_field = 10.0;
    zdplaskinSetConditions(&temperature, &ruduced_field);
}

}
