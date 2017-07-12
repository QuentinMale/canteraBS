/**
 * @file PlasmaKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAKINETICS_H
#define CT_PLASMAKINETICS_H

#include <Python.h>
#include "GasKinetics.h"
#include "RxnRates.h"
#include "Reaction.h"

using namespace std;

namespace Cantera
{
// forward declaration
class PlasmaRate;
class PlasmaReaction;

/**
 * Kinetics manager for plasma chemistry.
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
    PlasmaKinetics(thermo_t* thermo = 0);
    ~PlasmaKinetics(){
        Py_DECREF(m_processes);
        Py_DECREF(m_module);
        Py_DECREF(m_boltzmann);
        Py_DECREF(m_eedf);
    }
    virtual std::string kineticsType() const {
        return "Plasma";
    };

    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    virtual void updateROP();
    virtual void init();

    void update_EEDF();
    double getPlasmaReactionRate(size_t i);
    void PrintTotalRefCount();

protected:
    void addPlasmaReaction(PlasmaReaction& r);
    void modifyPlasmaReaction(size_t i, PlasmaReaction& r);
    Rate1<PlasmaRate> m_plasma_rates;
    //! Reaction index of each plasma reaction
    vector<size_t> m_plasmaIndex;
    vector<string> m_equations;
    vector<size_t> m_list;
    vector<double> m_x;

    PyObject *m_processes;
    PyObject *m_boltzmann;
    PyObject *m_eedf;
    PyObject *m_module;
};

class PlasmaRate //not useful
{
public:
    //! return the rate coefficient type.
    static int type() {
        return PLASMA_REACTION_RATECOEFF_TYPE;
    }

    PlasmaRate();
};

class PlasmaReaction : public Reaction
{
public:
    PlasmaReaction();
    PlasmaReaction(const Composition& reactants, const Composition& products,
                   const PlasmaRate& rate);
    PlasmaRate rate;
};

}

#endif