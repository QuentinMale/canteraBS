/**
 * @file PlasmaKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAKINETICS_H
#define CT_PLASMAKINETICS_H

#include "GasKinetics.h"
#include "cantera/plasma/PlasmaPhase.h"

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
    PlasmaKinetics(thermo_t* thermo=0);

    virtual std::string kineticsType() const {
        return "Plasma";
    }

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    virtual void init();
    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    virtual void invalidateCache();
    //@}

    //! Update rates of the plasma reaction
    virtual void update_rates_T();

protected:
    void addElectronTemperatureReaction(ElectronTemperatureReaction& r);
    void addPlasmaReaction(PlasmaReaction& r);
    void modifyElectronTemperatureReaction(size_t i, ElectronTemperatureReaction& r);

    //! Gas temperature.
    double m_temp_gas;

    //! Electron temperature for electron-temperature reactions
    double m_temp_e;

    //! Reaction index of each plasma reaction
    std::vector<size_t> m_plasmaIndex;

    //! Process index of each plasma reaction;
    std::vector<size_t> m_plasmaProcessIndex;

    //! Pointer to the plasma phase
    PlasmaPhase* m_plasma;
};

}

#endif
