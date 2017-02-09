//! @file IonFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_IONFLOW_H
#define CT_IONFLOW_H

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/oneD/StFlow.h"

namespace Cantera
{
/**
 * A class for ion flow.
 * @ingroup onedim
 */
class IonFlow : public FreeFlame
{
public:
    IonFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);
    //! set the solving phase
    virtual void setSolvingPhase(const size_t phase);
    //! set electric voltage at inlet and outlet
    virtual void setElectricVoltage(const double v1, const double v2);

    virtual void eval(size_t jg, doublereal* xg,
              doublereal* rg, integer* diagg, doublereal rdt);

    virtual void resize(size_t components, size_t points);

    virtual void _finalize(const doublereal* x);
    //! set to solve Poisson's equation on a point
    void solvePoissonEqn(size_t j=npos);
    //! set to fix voltage on a point
    void fixElectricPotential(size_t j=npos);
    bool doPoisson(size_t j) {
        return m_do_poisson[j];
    }
    //! set to solve velocity on a point
    void solveVelocity(size_t j=npos);
    //! set to fix velocity on a point
    void fixVelocity(size_t j=npos);
    bool doVelocity(size_t j) {
        return m_do_velocity[j];
    }

protected:
    virtual void updateTransport(doublereal* x, size_t j0, size_t j1);
    virtual void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    //! evaluate the residual for Poisson's equation 
    virtual void evalPoisson(size_t j, doublereal* x, doublereal* r, integer* diag, doublereal rdt);
    //! Solving phase one: the fluxes of charged species are turned off
    virtual void phaseOneDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    //! Solving phase two: the Prager's ambipolar-diffusion model is used       
    virtual void phaseTwoDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    //! Solving phase three: the Poisson's equation is added coupled by the electrical drift 
    virtual void phaseThreeDiffFluxes(const doublereal* x, size_t j0, size_t j1);
    //! flag for solving poisson's equation or not
    std::vector<bool> m_do_poisson;
    //! flag for solving the velocity or not
    std::vector<bool> m_do_velocity;

    // !electrical properties
    vector_int m_speciesCharge;

    // !index of species with charges
    std::vector<size_t> m_kCharge;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! mobility
    vector_fp m_mobi;

    //! solving phase
    int m_phase;

    //! The voltage 
    double m_inletVoltage;
    double m_outletVoltage;

    //! index of electron
    size_t m_kElectron;

    //! fixed electric potential value
    vector_fp m_fixedElecPoten;

    //! fixed velocity value
    vector_fp m_fixedVelocity;

    //! The fixed electric potential value at point j
    doublereal phi_fixed(size_t j) const {
        return m_fixedElecPoten[j];
    }

    //! The fixed velocity value at point j
    doublereal u_fixed(size_t j) const {
        return m_fixedVelocity[j];
    }

    //! electric potential
    doublereal phi(const doublereal* x, size_t j) const {
        return x[index(c_offset_P, j)];
    }  

    //! electric field
    doublereal E(const doublereal* x, size_t j) const {
        return -(phi(x,j+1)-phi(x,j))/(z(j+1)-z(j));
    }

    doublereal dEdz(const doublereal* x, size_t j) const {
        return 2*(E(x,j)-E(x,j-1))/(z(j+1)-z(j-1));
    }

    //! number density
    doublereal ND(const doublereal* x, size_t k, size_t j) const {
        return Avogadro * m_rho[j] * Y(x,k,j) / m_wt[k];
    }
};

}

#endif
