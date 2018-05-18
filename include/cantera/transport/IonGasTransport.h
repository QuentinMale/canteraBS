/**
 * @file GasTransport.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ION_GAS_TRANSPORT_H
#define CT_ION_GAS_TRANSPORT_H

#include "MixTransport.h"

namespace Cantera
{

class IonGasTransport : public MixTransport
{
public:
    //! Default constructor.
    IonGasTransport();

    virtual std::string transportType() const {
        return "Ion";
    }

    virtual void init(thermo_t* thermo, int mode, int log_level);
    virtual double viscosity();
    virtual double thermalConductivity();
    //! The ion mobilities are calculated by Blanc's law
    /*! Chiflikian, R. V. "The analog of Blanc’s law for drift velocities
     *  of electrons in gas mixtures in weakly ionized plasma."
     *  Physics of Plasmas 2.10 (1995): 3902-3909.
     */
    virtual void getMobilities(double* const mobi);
    //! The mixture transport for ionized gas
    //! The binary transport between two charged species is neglected.
    virtual void getMixDiffCoeffs(double* const d);

protected:
    //! setup parameters for n64 model
    void setupN64();
    virtual void fitDiffCoeffs(MMCollisionInt& integrals);
    double omega11_n64(const double tstar, const double gamma);

    //! electrical properties
    vector_int m_speciesCharge;

    //! index of species with charges
    std::vector<size_t> m_kIon;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! index of electron
    size_t m_kElectron;

    DenseMatrix m_gamma;
};

}

#endif
