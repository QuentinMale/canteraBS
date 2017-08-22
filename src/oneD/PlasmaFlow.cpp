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
    m_elec_field(2.5e6),
    m_elec_frequency(1e10),
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
            cout << speciesName << " " << k << endl;
        }
    }
    m_elec_power.resize(m_points, 0.0);
    m_ND.resize(m_points, 0.0);
    m_import_electron_transport = true;
}

void PlasmaFlow::resize(size_t components, size_t points) {
    IonFlow::resize(components, points);
    m_elec_power.resize(m_points, 0.0);
    m_ND.resize(m_points, 0.0);
}

void PlasmaFlow::evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax)
{
    IonFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);

    if (m_do_plasma) {
        for (size_t j = jmin; j <= jmax; j++) {
            if (j != 0 && grid(j) <= m_zfixed) {
            // if (grid(j) == m_zfixed) {
                // set gas
                setGas(x,j);
                for (size_t k = 0; k < m_nsp; k++) {
                    m_ND[k] = ND(x,k,j);
                }
                m_ND_t = ND_t(j);
                updateEEDF(x, j);
                double* wdot_plasma = NULL;
                zdplaskinGetPlasmaSource(&wdot_plasma);
                for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
                    size_t k = m_plasmaSpeciesIndex[i];
                    // multiply by the multiplier
                    // wdot_plasma[i] *= m_plasma_multiplier;
                    rsd[index(c_offset_Y + k, j)] += m_wt[k] * wdot_plasma[i] / m_rho[j];
                }

                // update electron power
                const double number_density = ND_t(j);
                m_elec_power[j] = zdplaskinGetElecPower(&number_density);
                rsd[index(c_offset_T, j)] += m_elec_power[j]
                                             * m_elec_num_density
                                             / (m_rho[j] * m_cp[j])
                                             * m_plasma_multiplier;
            }
        }
    }
}

void PlasmaFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x, j0, j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_kElectron != npos) {
            if (m_do_plasma) {
                // set number density to mid point
                for (size_t k = 0; k < m_nsp; k++) {
                    m_ND[k] = 0.5*(ND(x,k,j+1) + ND(x,k,j));
                }
                m_ND_t = 0.5 * (ND_t(j+1) + ND_t(j));
                updateEEDF(x,j);
                m_mobility[m_kElectron+m_nsp*j] = zdplaskinGetElecMobility(&m_ND_t);
                m_diff[m_kElectron+m_nsp*j] = zdplaskinGetElecDiffCoeff();
            } else {
                m_mobility[m_kElectron+m_nsp*j] = 0.4;
                m_diff[m_kElectron+m_nsp*j] = 0.4 * Boltzmann / ElectronCharge;
                m_diff[m_kElectron+m_nsp*j] *= m_thermo->temperature();
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

void PlasmaFlow::setElecField(double elec_field)
{
    m_elec_field = elec_field;
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
    for (size_t k : m_plasmaSpeciesIndex) {
        const char* species = m_thermo->speciesName(k).c_str();
        zdplaskinSetDensity(species, &m_ND[k]);
    }

    // set electron density
    if (m_kElectron != npos) {
        zdplaskinSetDensity("E", &m_ND[m_kElectron]);
    } else {
        zdplaskinSetDensity("E", &m_elec_num_density);
    }
    const double temperature = m_thermo->temperature();
    zdplaskinSetConditions(&temperature, &m_elec_field, &m_elec_frequency, &m_ND_t);
}

}
