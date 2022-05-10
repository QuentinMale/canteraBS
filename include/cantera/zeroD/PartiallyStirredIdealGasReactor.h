//! @file PartiallyStirredIdealGasReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PARTIALLYSTIRREDIDEALGASREACTOR_H
#define CT_PARTIALLYSTIRREDIDEALGASREACTOR_H

#include "IdealGasReactor.h"

namespace Cantera
{

/**
 * Class PartiallyStirredIdealGasReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 */
class PartiallyStirredIdealGasReactor : public IdealGasReactor
{
public:
    PartiallyStirredIdealGasReactor() {}

    virtual std::string type() const {
        return "PartiallyStirredIdealGasReactor";
    }

    void setMixingFactor(double mf) {
        m_mixingFactor = mf;
    }

    virtual void eval(double t, double* LHS, double* RHS);

protected:
    //! Mixing factor
    double m_mixingFactor;
};

}

#endif
