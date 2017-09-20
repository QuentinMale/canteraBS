/**
 * @file GasTransport.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ION_GAS_TRANSPORT_H
#define CT_ION_GAS_TRANSPORT_H

#include "GasTransport.h"

namespace Cantera
{

class IonGasTransport : public GasTransport
{
public:

protected:
    IonGasTransport(ThermoPhase* thermo=0);
    virtual void updateDiff_T();
    // get diffusion coefficient for charged-charged interaction
    double getCoulombDiffusion(const size_t i, const size_t j);

    //! electrical properties
    vector_int m_speciesCharge;

    //! index of species with charges
    std::vector<size_t> m_kCharge;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! index of electron
    size_t m_kElectron;

    //! Current value of the mean molecular weight
    double m_mmw;

    //! Current value of the density
    double m_rho;
};

}

#endif