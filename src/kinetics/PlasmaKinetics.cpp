/**
 *  @file PlasmaKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.
#include <Python.h>
#include "cantera/kinetics/PlasmaKinetics.h"

using namespace std;

namespace Cantera
{
PlasmaKinetics::PlasmaKinetics(thermo_t* thermo) :
    GasKinetics(thermo)
{
    Py_Initialize();
    PyList_Append(PySys_GetObject((char*)"path"),
        PyString_FromString("/home/bang/cantera/build/python2/cantera"));
    PyList_Append(PySys_GetObject((char*)"path"),
        PyString_FromString("/home/chen1671/cantera/build/python2/cantera"));
}

void PlasmaKinetics::calculateEEDF()
{
    //gas composition
    const vector<string> name = thermo().speciesNames();
    vector_fp x(thermo().nSpecies(), 0.0);
    thermo().getMoleFractions(&x[0]);

    // python object
    PyObject *pModule, *pFunc;
    PyObject *pArgs;
    PyObject *gas_species, *gas_molefraction;

    gas_species = PyList_New(0);
    gas_molefraction = PyList_New(0);
    for (size_t i = 0; i < thermo().nSpecies(); i++) {
        PyList_Append(gas_species,PyString_FromString(name[i].c_str()));
        PyList_Append(gas_molefraction,PyFloat_FromDouble(x[i]));
    }

    const char *fileName = "eedf";
    const char *funcName = "eedf";

    pModule = PyImport_Import(PyString_FromString(fileName));
    pFunc = PyObject_GetAttrString(pModule, funcName);
    pArgs = PyTuple_New(3);
    PyTuple_SetItem(pArgs, 0, gas_species);
    PyTuple_SetItem(pArgs, 1, gas_molefraction);
    PyTuple_SetItem(pArgs, 2, PyFloat_FromDouble(thermo().temperature()));
    PyObject_CallObject(pFunc, pArgs);

    Py_DECREF(pModule);
    Py_DECREF(pFunc);
    Py_DECREF(pArgs);
    Py_DECREF(gas_species);
    Py_DECREF(gas_molefraction);
}

double PlasmaKinetics::getPlasmaReactionRate(string equation)
{
    cout << equation << endl;
    return 87;
}

void PlasmaKinetics::updateROP()
{
    GasKinetics::updateROP();
    calculateEEDF();
    vector_fp pr(m_plasma_rates.nReactions(),0.0);
    for (size_t i = 0; i < m_plasma_rates.nReactions(); i++) {
        pr[i] = getPlasmaReactionRate(m_equations[i]);

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
