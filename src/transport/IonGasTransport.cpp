//! @file IonGasTransport.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/IonGasTransport.h"
#include "cantera/zdplaskin.h"

namespace Cantera
{
IonGasTransport::IonGasTransport(ThermoPhase* thermo) :
    GasTransport(thermo),
    m_mmw(0.0),
    m_kElectron(npos)
{
    m_thermo = thermo;
    m_nsp = m_thermo->nSpecies();
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            m_kCharge.push_back(k);
        } else {
            m_kNeutral.push_back(k);
        }
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }
}

void IonGasTransport::updateDiff_T()
{
    GasTransport::updateDiff_T();
    update_C();
    // evaluate binary diffusion coefficients for charge-charge interaction
    // update mean molecular weight
    m_mmw = m_thermo->meanMolecularWeight();
    // update density
    m_rho = m_thermo->density();
    for (size_t i : m_kCharge) {
        for (size_t j : m_kCharge) {
            m_bdiff(i,j) = getCoulombDiffusion(i,j);
            //cout << "charge-charge diff: " << m_bdiff(i,j) << endl;
        }
    }
    m_bindiff_ok = true;
}

double IonGasTransport::getCoulombDiffusion(const size_t i, const size_t j)
{
    const double reducedMass = m_mw[i] * m_mw[j] / (m_mw[i] + m_mw[j]);
    const double number_density = Avogadro * m_rho / m_mmw;
    double sum = 0.0;
    for (size_t k : m_kCharge) {
        sum += number_density * m_molefracs[k] * m_speciesCharge[k] * m_speciesCharge[k];
    }
    // the collision diameter is equal to debye length
    double lambda_D = sqrt(epsilon_0 * Boltzmann * m_temp /
                           (ElectronCharge * ElectronCharge * sum));
    // The well depth
    double epsilon = abs(m_speciesCharge[i] * m_speciesCharge[j]) * ElectronCharge * ElectronCharge
                     / (4 * Pi * epsilon_0 * lambda_D);

    // The collision integral is calculated using the fitting curve in references:
    // Han, Jie, et al. "Numerical modelling of ion transport in flames."
    // Combustion Theory and Modelling 19.6 (2015): 744-772.
    double tstar = Boltzmann * m_temp / epsilon;
    // reduced collision integral
    double om11 = (0.5*log(tstar) - 0.14) / (tstar * tstar);
    // cross section
    double sigma = Pi * lambda_D * lambda_D * om11;

    double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi * Avogadro * Boltzmann * m_temp / reducedMass)
                       / (sigma * number_density);

    return diffcoeff;
}

}
