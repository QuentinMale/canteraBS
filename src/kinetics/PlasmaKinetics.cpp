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

void PlasmaKinetics::addPlasmaReaction(PlasmaReaction& r)
{
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

void PlasmaKinetics::init()
{
    GasKinetics::init();
    m_plasma = dynamic_cast<PlasmaPhase*>(&thermo());
}

}
