//! @file PlasmaReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAREACTOR_H
#define CT_PLASMAREACTOR_H

#include "IdealGasReactor.h"
#include "cantera/thermo/PlasmaPhase.h"

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

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void updateState(double* y) override;

    void eval(double t, double* LHS, double* RHS) override;

    //! Set/Get discharge volume
    void setDisVol(double dis_vol) {
        m_dis_vol = dis_vol;
    }
    double disVol() const {
        return m_dis_vol;
    }

    //! Get discharge volumetric power
    //CQM may not be up to date
    double disVPower() const{
        return m_disVPower;
    }

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

    void compute_disVPower();

protected:
    void setThermo(ThermoPhase& thermo) override;

    double m_dis_vol; //!< Discharge volume

    double m_disVPower; //!< Volumetric discharge power

    size_t m_nspevib = 1; //!< Number of species with vibrational excitation

    PlasmaPhase* m_plasma = nullptr;
};
}

#endif
