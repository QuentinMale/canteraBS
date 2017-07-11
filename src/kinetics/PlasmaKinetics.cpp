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
    // Dangerous! don't touch anything!
    Py_Initialize();
    PyList_Append(PySys_GetObject((char*)"path"),
        PyString_FromString("."));

    m_module = PyImport_Import(Py_BuildValue("s","eedf"));
    if (!m_module)
    {
        cout << "Error importing bolos interface" << endl;
    }
    PyObject *pFunc = PyObject_GetAttrString(m_module, "initialize");

    m_processes = PyObject_CallObject(pFunc, NULL);
    Py_DECREF(pFunc);
}

void PlasmaKinetics::update_EEDF()
{
    // Dangerous! don't touch anything!
    //gas composition
    const vector<string> name = thermo().speciesNames();
    vector_fp x(thermo().nSpecies(), 0.0);
    thermo().getMoleFractions(&x[0]);

    PyObject *pGasSpecies = PyList_New(0);
    PyObject *pGasMoleFraction = PyList_New(0);
    for (size_t i = 0; i < thermo().nSpecies(); i++) {
        PyObject *pName = Py_BuildValue("s",name[i].c_str());
        PyList_Append(pGasSpecies, pName);
        Py_DECREF(pName);
        PyObject *pX = Py_BuildValue("d",x[i]);
        PyList_Append(pGasMoleFraction, pX);
        Py_DECREF(pX);
    }
    PyObject *pTemp = Py_BuildValue("d",thermo().temperature());
    PyObject *pFunc = PyObject_GetAttrString(m_module, "eedf");
    PyObject *ptuple = PyObject_CallFunctionObjArgs(pFunc, 
                                                    m_processes,
                                                    pGasSpecies, 
                                                    pGasMoleFraction, 
                                                    pTemp,
                                                    NULL);
    Py_DECREF(pFunc);
    Py_DECREF(pGasSpecies);
    Py_DECREF(pGasMoleFraction);
    Py_DECREF(pTemp);
    m_eedf = Py_BuildValue("O", PyTuple_GetItem(ptuple, 0));
    m_boltzmann = Py_BuildValue("O", PyTuple_GetItem(ptuple, 1));
    Py_DECREF(ptuple);
}

double PlasmaKinetics::getPlasmaReactionRate(size_t i)
{
    // Dangerous! don't touch anything!
    PyObject *pEquation = Py_BuildValue("s",m_equations[i].c_str());
    PyObject *pFunc = PyObject_GetAttrString(m_module, "getReactionRate");
    PyObject *pK = PyObject_CallFunctionObjArgs(pFunc, 
                                                m_eedf, 
                                                m_boltzmann,
                                                pEquation,
                                                NULL);
    Py_DECREF(pFunc);
    Py_DECREF(pEquation);
    double k = PyFloat_AsDouble(pK);
    Py_DECREF(pK);
    return k;
}
/*
void PlasmaKinetics::PrintTotalRefCount()
{
#ifdef Py_REF_DEBUG
    PyObject* refCount = PyObject_CallObject(PySys_GetObject((char*)"gettotalrefcount"), NULL);
    std::clog << "total refcount = " << PyInt_AsSsize_t(refCount) << std::endl;
    Py_DECREF(refCount);
#endif
}
*/
void PlasmaKinetics::updateROP()
{
    GasKinetics::updateROP();
    update_EEDF();
    vector_fp pr(m_plasma_rates.nReactions(),0.0);
    for (size_t i = 0; i < m_plasma_rates.nReactions(); i++) {
        pr[i] = getPlasmaReactionRate(i);

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
