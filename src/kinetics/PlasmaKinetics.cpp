/**
 *  @file PlasmaKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.
#include "cantera/kinetics/PlasmaKinetics.h"

using namespace std;

namespace Cantera
{
PlasmaKinetics::PlasmaKinetics(thermo_t* thermo) :
    GasKinetics(thermo),
    m_do_plasma(true)
{
}

void PlasmaKinetics::init()
{
    GasKinetics::init();
    vector<string> list;
    list.push_back("N2");
    list.push_back("N");
    list.push_back("O2");
    list.push_back("O");
    list.push_back("H2");
    list.push_back("H");
    list.push_back("CO2");
    list.push_back("CO");
    list.push_back("H2O");
    list.push_back("CH4");
    // check valid index 
    for (size_t i = 0; i < list.size(); i++) {
        size_t k = kineticsSpeciesIndex(list[i]);
        if (k != npos) {
            m_list.push_back(k);
        }
    }

    m_x.resize(thermo().nSpecies(), 0.0);
}

void PlasmaKinetics::update_rates_C()
{
    GasKinetics::update_rates_C();
    // update plasma reaction rates
    if (m_do_plasma) {
        // need to update eedf somewhere
        vector_fp pr(m_plasma_rates.nReactions(),0.0);
        for (size_t i = 0; i < m_plasma_rates.nReactions(); i++) {
            // need to update eedf somewhere
            //pr[i] = m_plasma_rate_coeffs[i];
            pr[i] = 0.0;

            AssertFinite(pr[i], "PlasmaKinetics::updateROP",
                         "pr[{}] is not finite.", i);
        }
        scatter_copy(pr.begin(), pr.begin() + m_plasma_rates.nReactions(),
                     m_rfn.begin(), m_plasmaIndex.begin());
    }
}

void PlasmaKinetics::solvePlasmaRates(bool doPlasma)
{
    m_do_plasma = doPlasma;
}

bool PlasmaKinetics::addReaction(shared_ptr<Reaction> r)
{
    if (r->reaction_type == PLASMA_RXN) {
        // operations common to all reaction types
        bool added = BulkKinetics::addReaction(r);
        if (!added) {
            return false;
        }
        addPlasmaReaction(dynamic_cast<PlasmaReaction&>(*r));
        return true;
    } else {
        bool added = GasKinetics::addReaction(r);
        return added;
    }
}

void PlasmaKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    if (rNew->reaction_type == PLASMA_RXN) {
        modifyPlasmaReaction(i, dynamic_cast<PlasmaReaction&>(*rNew));
    } else {
        GasKinetics::modifyReaction(i, rNew);
    }
}

void PlasmaKinetics::addPlasmaReaction(PlasmaReaction& r)
{
    m_equations.push_back(r.equation());
    m_plasmaIndex.push_back(nReactions()-1);
    m_plasma_rates.install(nReactions()-1, r.rate);
}

void PlasmaKinetics::modifyPlasmaReaction(size_t i, PlasmaReaction& r)
{
    m_plasma_rates.replace(i, r.rate);
}

// class PlasmaRate
PlasmaRate::PlasmaRate()
    // composition
{
}

// class PlasmaReaction
PlasmaReaction::PlasmaReaction() :
    Reaction(PLASMA_RXN)
{
}

PlasmaReaction::PlasmaReaction(const Composition& reactants_,
                               const Composition& products_,
                               const PlasmaRate& rate_) :
    Reaction(PLASMA_RXN, reactants_, products_),
    rate(rate_)
{
}

}
