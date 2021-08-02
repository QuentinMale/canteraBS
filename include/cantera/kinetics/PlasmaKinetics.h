/**
 * @file PlasmaKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAKINETICS_H
#define CT_PLASMAKINETICS_H

#include "GasKinetics.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary plasma-phase chemistry.
 * @ingroup kinetics
 */
class PlasmaKinetics : public GasKinetics
{
public:
    //! @name Constructors and General Information
    //! @{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    PlasmaKinetics(ThermoPhase* thermo=0);

    virtual std::string kineticsType() const {
        return "Plasma";
    }

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    // virtual void init();
    virtual void invalidateCache();
    //! @}

    //! Update rates of the plasma reaction
    virtual void update_rates_T();

protected:
    //! Add a specific type of reaction
    virtual bool addSpecificReaction(shared_ptr<Reaction> r);
    //! Modify a specific type of reaction
    virtual bool modifySpecificReaction(size_t i, shared_ptr<Reaction> r);
    //! Add a electron-temperature type of reaction
    void addETempReaction(ETempReaction& r);
    //! Modify a electron-temperature type of reaction
    void modifyETempReaction(size_t i, ETempReaction& r);

    //! Rate expressions for electron temperature reaction
    Rate1<ElectronTemperature> m_electron_temp_rates;

    //! Gas temperature.
    double m_temp_gas;

    //! Electron temperature.
    double m_temp_e;
};

}

#endif
