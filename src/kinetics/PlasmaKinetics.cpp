/**
 *  @file PlasmaKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PlasmaKinetics.h"
#include <Python.h>

using namespace std;

namespace Cantera
{
PlasmaKinetics::PlasmaKinetics(thermo_t* thermo) :
    GasKinetics(thermo)
{
    Py_Initialize();
    PyRun_SimpleString("print('This is python interpreter in PlasmaKinetics')\n");
    //Py_Finalize();
}

void PlasmaKinetics::updateROP()
{
    GasKinetics::updateROP();
    vector_fp pr(m_plasma_rates.nReactions(),0.0);
    for (size_t i = 0; i < m_plasma_rates.nReactions(); i++) {
        pr[i] = 87;
        AssertFinite(pr[i], "PlasmaKinetics::updateROP",
                     "pr[{}] is not finite.", i);
    }
    scatter_copy(pr.begin(), pr.begin() + m_plasma_rates.nReactions(),
                 m_ropf.begin(), m_plasmaIndex.begin());
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
