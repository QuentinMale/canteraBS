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

    virtual void resize(size_t components, size_t points);
    void updateEEDF(double* x, size_t j);
    void solvePlasma();
    bool doPlasma() {
        return m_do_plasma;
    }
    void setElecField(double reduced_field);
    void setElecNumDensity(double num_density);
    void setPlasmaSourceMultiplier(double multiplier);

protected:
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);
    virtual void updateTransport(double* x, size_t j0, size_t j1);
    vector<size_t> m_collisionSpeciesIndex;
    vector<size_t> m_plasmaSpeciesIndex;
    bool m_do_plasma;
    double m_elec_num_density;
    double m_reduced_field;
    double m_plasma_multiplier;
};

}

#endif
