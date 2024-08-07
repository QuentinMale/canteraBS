//! @file PlasmaReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAREACTOR_H
#define CT_PLASMAREACTOR_H

#include "IdealGasReactor.h"

namespace Cantera
{

/**
 * Class PlasmaReactor is a class for stirred reactors that ...
 *
 *
 * @ingroup reactorGroup
 */
class PlasmaReactor : public IdealGasReactor
{
public:
    using IdealGasReactor::IdealGasReactor; // inherit constructors

    string type() const override {
        return "PlasmaReactor";
    }

    //! Set/Get discharge volume
    void setDisVol(doublereal dis_vol) {
        m_dis_vol = dis_vol;
    }
    double disVol() {
        return m_dis_vol;
    }

protected:
    void setThermo(ThermoPhase& thermo) override;

    double m_dis_vol; //!< Discharge volume

};
}

#endif
