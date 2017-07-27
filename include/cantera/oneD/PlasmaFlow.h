//! @file PlasmaFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAFLOW_H
#define CT_PLASMAFLOW_H

#include "cantera/oneD/IonFlow.h"

using namespace std;

namespace Cantera
{
/**
 * This class models plasma in a flame. There are three
 * stages of the simulation.
 * @ingroup onedim
 */
class PlasmaFlow : public IonFlow
{
public:
    PlasmaFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    virtual void eval(size_t jg, double* xg,
              double* rg, integer* diagg, double rdt);
    virtual void resize(size_t components, size_t points);
    void updateEEDF();

protected:
    virtual void updateTransport(double* x, size_t j0, size_t j1);
    vector<size_t> m_collisionSpeciesIndex;
    vector<size_t> m_plasmaSpeciesIndex;

    virtual void getWdot(doublereal* x, size_t j);
    double* m_wdot_plasma = NULL;
    // ZDPlasKin wrappe

};

}

#endif
