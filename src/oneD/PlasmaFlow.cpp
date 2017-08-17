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
    m_do_plasma(false),
    m_elec_num_density(1e17),
    m_reduced_field(100),
    m_plasma_multiplier(1.0)
{
    vector<string> list;
    list.push_back("N2");
    list.push_back("O2");
    list.push_back("H2");
    list.push_back("H2O");
    list.push_back("CO2");
    list.push_back("CO");
    list.push_back("CH4");
    list.push_back("E");
    // check valid index 
    for (size_t i = 0; i < list.size(); i++) {
        size_t k = m_thermo->speciesIndex(list[i]);
        if (k != npos) {
            m_collisionSpeciesIndex.push_back(k);
        }
    }

    for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
        char cstring[20];
        zdplaskinGetSpeciesName(cstring, &i);
        string speciesName(cstring);
        size_t k = m_thermo->speciesIndex(speciesName);
        if (k != npos) {
            m_plasmaSpeciesIndex.push_back(k);
        }
    }
    m_elec_power.resize(m_points, 0.0);
}

void PlasmaFlow::resize(size_t components, size_t points) {
    IonFlow::resize(components, points);
    m_elec_power.resize(m_points, 0.0);
}

void PlasmaFlow::evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax)
{
    IonFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);
    size_t j0 = max<size_t>(jmin, 1);
    size_t j1 = min(jmax, m_points-2);

    if (m_do_plasma) {
        m_do_plasma = false;
        for (size_t j = j0; j <= j1; j++) {
            updateEEDF(x, j);
            // double* wdot_plasma = NULL;
            // zdplaskinGetPlasmaSource(&wdot_plasma);
            // for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
            //     size_t k = m_plasmaSpeciesIndex[i];
            //     // multiply by the multiplier
            //     wdot_plasma[i] *= m_plasma_multiplier;
            //     rsd[index(c_offset_Y + k, j)] += m_wt[k] * wdot_plasma[i] / m_rho[j];
            // }

            // update electron power
            const double number_density = ND_t(j);
            m_elec_power[j] = zdplaskinGetElecPower(&number_density);
        }
    }

    // add electron power to the system
    for (size_t j = j0; j <= j1; j++) {
        rsd[index(c_offset_T, j)] += m_elec_power[j]
                                     * m_elec_num_density 
                                     / (m_rho[j] * m_cp[j])
                                     * m_plasma_multiplier;
    }
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
                const double number_density = ND_t(j);
                m_mobility[m_kElectron+m_nsp*j] = zdplaskinGetElecMobility(&number_density);
                m_diff[m_kElectron+m_nsp*j] = zdplaskinGetElecDiffCoeff();
            } else {
                m_mobility[m_kElectron+m_nsp*j] = 0.4;
                m_diff[m_kElectron+m_nsp*j] = 0.4*(Boltzmann * T(x,j)) / ElectronCharge;
            }
        }
    }
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

void PlasmaFlow::setElecField(double reduced_field)
{
    m_reduced_field = reduced_field;
}

void PlasmaFlow::setElecNumDensity(double num_density)
{
    m_elec_num_density = num_density;
}

void PlasmaFlow::setPlasmaSourceMultiplier(double multiplier)
{
    m_plasma_multiplier = multiplier; 
}

void PlasmaFlow::updateEEDF(double* x, size_t j)
{
    for (size_t k : m_collisionSpeciesIndex) {
        const char* species = m_thermo->speciesName(k).c_str();
        const double num_density = ND(x,k,j);
        zdplaskinSetDensity(species, &num_density);
    }

    // set electron density
    zdplaskinSetDensity("E", &m_elec_num_density);

    const double temperature = T(x,j);
    zdplaskinSetConditions(&temperature, &m_reduced_field);
}

}
