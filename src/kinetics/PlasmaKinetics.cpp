/**
 *  @file PlasmaKinetics.cpp Kinetics in plasmas
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PlasmaKinetics.h"
#include "cantera/thermo/ThermoPhase.h"

using namespace std;

namespace Cantera
{
PlasmaKinetics::PlasmaKinetics(ThermoPhase* thermo)
    : GasKinetics(thermo)
{
}

void PlasmaKinetics::update_rates_T()
{
    GasKinetics::update_rates_T();
    // update electron-temperture reactions
    double T = thermo().temperature();
    double Te = thermo().electronTemperature();
    double logTe = log(Te);
    if (T != m_temp_gas || Te != m_temp_e) {
        if (m_electron_temp_rates.nReactions()) {
            m_electron_temp_rates.updateETemp(T, Te, logTe, m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_temp_e = Te;
    m_temp_gas = T;
}

bool PlasmaKinetics::addSpecificReaction(shared_ptr<Reaction> r)
{
    if (GasKinetics::addSpecificReaction(r)){
        return true;
    } else if (r->type() == "electron-temperature") {
        addETempReaction(dynamic_cast<ETempReaction&>(*r));
    } else {
        return false;
    }
    return true;
}

void PlasmaKinetics::addETempReaction(ETempReaction& r)
{
    m_electron_temp_rates.install(nReactions()-1, r.rate);
}

bool PlasmaKinetics::modifySpecificReaction(size_t i, shared_ptr<Reaction> rNew)
{
    if (GasKinetics::modifySpecificReaction(i, rNew)) {
        return true;
    } else if (rNew->type() == "electron-temperature") {
        modifyETempReaction(i, dynamic_cast<ETempReaction&>(*rNew));
    } else {
        return false;
    }
    return true;
}

void PlasmaKinetics::modifyETempReaction(size_t i, ETempReaction& r)
{
    m_electron_temp_rates.replace(i, r.rate);
}

void PlasmaKinetics::invalidateCache()
{
    GasKinetics::invalidateCache();
    m_temp_e += 0.13579;
}

}
