/**
 *  @file PlasmaKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PlasmaKinetics.h"

using namespace std;

namespace Cantera
{
PlasmaKinetics::PlasmaKinetics(thermo_t* thermo)
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
        if (m_electron_temperature_rates.nReactions()) {
            m_electron_temperature_rates.update(T, Te, logTe, m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_temp_e = Te;
    m_temp_gas = T;

    // update plasma reactions
    auto i = 0;
    for (auto k : m_plasmaProcessIndex) {
        m_rfn[m_plasmaIndex[i]] = m_plasma->rateCoefficient(k);
        if (m_rkcn[m_plasmaIndex[i]] != 0.0) {
            m_rkcn[m_plasmaIndex[i]] = m_plasma->reverseRateCoefficient(k) /
                                      m_rfn[m_plasmaIndex[i]];
        }
        // convert to kmol base
        m_rfn[m_plasmaIndex[i]] *= Avogadro;
        i++;
    }
    if (i != 0) {
        m_ROP_ok = false;
    }
}

bool PlasmaKinetics::addReaction(shared_ptr<Reaction> r)
{
    // operations common to all reaction types
    bool added = GasKinetics::addReaction(r);
    if (added) {
        return true;
    }

    switch (r->reaction_type) {
    case ELECTRON_TEMPERATURE_RXN:
        addElectronTemperatureReaction(dynamic_cast<ElectronTemperatureReaction&>(*r));
        break;
    case PLASMA_RXN:
        addPlasmaReaction(dynamic_cast<PlasmaReaction&>(*r));
        break;
    default:
        throw CanteraError("PlasmaKinetics::addReaction",
            "Unknown reaction type specified: {}", r->reaction_type);
    }
    return true;
}

void PlasmaKinetics::addElectronTemperatureReaction(ElectronTemperatureReaction& r)
{
    m_electron_temperature_rates.install(nReactions()-1, r.rate);
}

void PlasmaKinetics::addPlasmaReaction(PlasmaReaction& r)
{
    if (m_plasma == 0) {
        throw CanteraError("PlasmaKinetics::addPlasmaReaction",
                           "Reaction {} requires a class derived from PlasmaPhase"
                           "as thermo. Ex. thermo: weakly-ionized-gas", r.equation());
    }
    m_plasmaIndex.push_back(nReactions()-1);

    // find plasma process index
    bool found = false;
    for (size_t k = 0; k < m_plasma->nElectronCrossSections(); k++) {
        if (r.process_kind == m_plasma->kind(k) &&
            r.process_target == m_plasma->target(k) &&
            r.process_product == m_plasma->product(k)) {
            m_plasmaProcessIndex.push_back(k);
            found = true;
            break;
        }
    }
    if (!found) {
        throw CanteraError("PlasmaKinetics::addPlasmaReaction",
                           "Cannot find corresponding electron "
                           "collision process for {}", r.equation());
    }
}

void PlasmaKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    GasKinetics::modifyReaction(i, rNew);
    switch (rNew->reaction_type) {
    case ELECTRON_TEMPERATURE_RXN:
        modifyElectronTemperatureReaction(i, dynamic_cast<ElectronTemperatureReaction&>(*rNew));
        break;
    case PLASMA_RXN:
        throw CanteraError("GasKinetics::modifyReaction",
            "{} reaction type cannot be modified", rNew->reaction_type);
        break;
    default:
        throw CanteraError("GasKinetics::modifyReaction",
            "Unknown reaction type specified: {}", rNew->reaction_type);
    }
    m_temp_e += 0.1234;
}

void PlasmaKinetics::modifyElectronTemperatureReaction(size_t i, ElectronTemperatureReaction& r)
{
    m_electron_temperature_rates.replace(i, r.rate);
}


void PlasmaKinetics::init()
{
    GasKinetics::init();
    m_plasma = dynamic_cast<PlasmaPhase*>(&thermo());
}

void PlasmaKinetics::invalidateCache()
{
    GasKinetics::invalidateCache();
    m_temp_e += 0.13579;
}

}
