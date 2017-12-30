//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/zdplaskin.h"

using namespace std;

namespace Cantera
{

IonFlow::IonFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    FreeFlame(ph, nsp, points),
    m_do_elec_heat(false),
    m_stage(1),
    m_inletVoltage(0.0),
    m_outletVoltage(0.0),
    m_kElectron(npos),
    m_elec_field(0.0),
    m_elec_frequency(0.0),
    m_plasma_multiplier(1.0),
    m_electron_multiplier(1.0)
{
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            m_kCharge.push_back(k);
        } else {
            m_kNeutral.push_back(k);
        }
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
        char cstring[20];
        zdplaskinGetSpeciesName(cstring, &i);
        string speciesName(cstring);
        size_t k = m_thermo->speciesIndex(speciesName);
        m_kPlasmaSpecies.push_back(k);
    }

    // collision list
    vector<string> collision_list;
    // collision_list.push_back("N2");
    collision_list.push_back("O2");
    collision_list.push_back("CH4");
    collision_list.push_back("CO2");
    collision_list.push_back("H2O");
    collision_list.push_back("H2");
    collision_list.push_back("CO");

    for (size_t i = 0; i < collision_list.size(); i++) {
        size_t k = m_thermo->speciesIndex(collision_list[i]);
        if (k != npos ) {
            m_kCollision.push_back(k);
        }
    }

    // no bound for electric potential
    setBounds(c_offset_P, -1.0e20, 1.0e20);

    m_refiner->setActive(c_offset_P, false);
    m_mobility.resize(m_nsp*m_points);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
    m_do_plasma.resize(m_points,false);
    m_electronPower.resize(m_points, 0.0);
    m_electronTemperature.resize(m_points, 0.0);
    m_electronMobility.resize(m_points, 0.0);
    m_electronDiff.resize(m_points, 0.0);
    m_Eambi.resize(m_points, 0.0);
    m_wdotPlasma.resize(zdplaskinNSpecies(),m_points,0.0);
}

void IonFlow::resize(size_t components, size_t points){
    StFlow::resize(components, points);
    m_mobility.resize(m_nsp*m_points);
    m_do_species.resize(m_nsp,true);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
    m_do_plasma.resize(m_points,false);
    m_fixedElecPoten.resize(m_points,0.0);
    m_fixedVelocity.resize(m_points);
    m_electronPower.resize(m_points, 0.0);
    m_electronTemperature.resize(m_points, 0.0);
    m_electronMobility.resize(m_points, 0.0);
    m_electronDiff.resize(m_points, 0.0);
    m_Eambi.resize(m_points, 0.0);
    m_wdotPlasma.resize(zdplaskinNSpecies(),m_points,0.0);
}

void IonFlow::updatePlasmaProperties(const double* x)
{
    for (size_t j = 0; j < m_points - 1; j++) {
        // update electron properties
        if (m_do_plasma[j]) {
            for (size_t k : m_kCollision) {
                // if (ND(x,k,j) < 0.0) {
                //     number_density = -ND(x,k,j);
                // }
                double number_density = ND(x,k,j);
                const char* species = m_thermo->speciesName(k).c_str();
                zdplaskinSetDensity(species, &number_density);
            }
            // set electron number density
            double number_density = ND(x,m_kElectron,j);
            zdplaskinSetDensity("E", &number_density);

            const double Tgas = T(x,j);
            double total_number_density = ND_t(j);
            if (m_elec_field == 0) {
                //set electron temperature if E = 0
                zdplaskinSetElecTemp(&Tgas);
            } else {
                zdplaskinSetElecField(&m_elec_field, &m_elec_frequency, &total_number_density);
            }
            //set gas temperature
            zdplaskinSetGasTemp(&Tgas);

            // get plasma properties
            double multi = m_electron_multiplier;
            m_electronTemperature[j] = zdplaskinGetElecTemp();
            m_electronMobility[j] = zdplaskinGetElecMobility(&total_number_density) * multi
                                    + 0.4 * (1.0 - multi);
            m_electronDiff[j] = zdplaskinGetElecDiffCoeff() * multi
                                + 0.4*(Boltzmann * T(x,j)) / ElectronCharge * (1.0 - multi);
            m_electronPower[j] = zdplaskinGetElecPowerElastic(&total_number_density) * multi;
            m_electronPower[j] += zdplaskinGetElecPowerInelastic(&total_number_density) * multi;

            double* wdot_plasma = NULL;
            zdplaskinGetPlasmaSource(&wdot_plasma);
            for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
                m_wdotPlasma(i,j) = wdot_plasma[i];
            }
        } else {
            m_electronTemperature[j] = T(x,j);
            m_electronMobility[j] = m_mobility[m_kElectron+m_nsp*j];
            m_electronDiff[j] = m_diff[m_kElectron+m_nsp*j];
            m_electronPower[j] = 0.0;
        }
    }
    m_electronTemperature[m_points-1] = m_electronTemperature[m_points-2];
    m_electronMobility[m_points-1] = m_electronMobility[m_points-2];
    m_electronDiff[m_points-1] = m_electronDiff[m_points-2];
    m_electronPower[m_points-1] = m_electronPower[m_points-2];
}

void IonFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x, j0, j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_kElectron != npos) {
            size_t k = m_kElectron;
            if (m_do_plasma[j]) {
                m_mobility[k+m_nsp*j] = 0.5*(m_electronMobility[j]+m_electronMobility[j+1]);
                m_diff[k+m_nsp*j] = 0.5*(m_electronDiff[j]+m_electronDiff[j+1]);
            } else {
                m_mobility[k+m_nsp*j] = 0.4;
                m_diff[k+m_nsp*j] = 0.4*(Boltzmann * T(x,j)) / ElectronCharge;
            }
        }
    }
}

void IonFlow::evalResidual(double* x, double* rsd, int* diag,
                           double rdt, size_t jmin, size_t jmax)
{
    StFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);
    if (m_stage == 3) {
        // update plasma properties
        updatePlasmaProperties(x);
        for (size_t j = jmin; j <= jmax; j++) {
            if (j == 0) {
                // force phi will result bad electron profile
                //rsd[index(c_offset_P, j)] = phi(x,j) - m_inletVoltage;
                size_t k = m_kElectron;
                rsd[index(c_offset_Y + k, 0)] = Y(x,k,0);
                rsd[index(c_offset_P, j)] = E(x,j);
                diag[index(c_offset_P, j)] = 0;
            } else if (j == m_points - 1) {
                rsd[index(c_offset_P, j)] = m_outletVoltage - phi(x,j);
                diag[index(c_offset_P, j)] = 0;
            } else {
                //-----------------------------------------------
                //    Poisson's equation
                //
                //    dE/dz = e/eps_0 * sum(q_k*n_k)
                //
                //    E = -dV/dz
                //-----------------------------------------------
                rsd[index(c_offset_P, j)] = dEdz(x,j) - rho_e(x,j) / epsilon_0;
                diag[index(c_offset_P, j)] = 0;

                // This method is used when you disable energy equation
                // but still maintain the velocity profile
                if (!m_do_velocity[j]) {
                    rsd[index(c_offset_U, j)] = u(x,j) - u_fixed(j);
                    diag[index(c_offset_U, j)] = 0;
                }
                if (m_do_plasma[j]) {
                    for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
                        size_t k = m_kPlasmaSpecies[i];
                        if (k != npos) {
                            // multiply by the multiplier
                            rsd[index(c_offset_Y + k, j)] += m_wt[k] * m_wdotPlasma(i,j) / m_rho[j] *
                                                             m_plasma_multiplier;
                        }
                    }

                    // update electron power
                    if (m_do_elec_heat) {
                        if (m_do_energy[j]) {
                            size_t k = m_kElectron;
                            double ND_e = (ND(x,k,j) > 0.0 ? ND(x,k,j) : 0.0);
                            rsd[index(c_offset_T, j)] += m_electronPower[j]
                                                         * ND_e / (m_rho[j] * m_cp[j])
                                                         * m_electron_multiplier;
                        }
                    }
                }
            }
        }
    }
}

void IonFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_stage == 1) {
        frozenIonMethod(x,j0,j1);
    }
    if (m_stage == 2) {
        chargeNeutralityModel(x,j0,j1);
    }
    if (m_stage == 3) {
        poissonEqnMethod(x,j0,j1);
    }
}

void IonFlow::frozenIonMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);
        double sum = 0.0;
        for (size_t k : m_kNeutral) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j);
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += sum*Y(x,k,j);
        }

        // flux for ions
        // Set flux to zero to prevent some fast charged species (e.g. electron)
        // to run away
        for (size_t k : m_kCharge) {
            m_flux(k,j) = 0;
        }
    }
}

void IonFlow::chargeNeutralityModel(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);

        // mixture-average diffusion
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        double sum_flux = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum_flux -= m_flux(k,j); // total net flux
        }
        double sum_ion = 0.0;
        for (size_t k : m_kCharge) {
            sum_ion += Y(x,k,j);
        }
        // The portion of correction for ions is taken off
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += Y(x,k,j) / (1-sum_ion) * sum_flux;
        }

        // ambipolar diffusion
        double sum_chargeFlux = 0.0;
        for (size_t k : m_kCharge) {
            double q_k = m_speciesCharge[k] * ElectronCharge;
            sum_chargeFlux += q_k * Avogadro / m_wt[k] * m_flux(k,j);
        }
        m_Eambi[j] = -sum_chargeFlux / sigma(x,j);
        for (size_t k : m_kCharge) {
            double Xav = 0.5 * (X(x,k,j+1) + X(x,k,j));
            double drift = s_k(k) * m_mobility[k+m_nsp*j] * m_Eambi[j] 
                           * rho * m_wt[k] / wtm * Xav;
            m_flux(k,j) += drift;
        }
    }
}

void IonFlow::poissonEqnMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);

        // mixture-average diffusion
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        double sum_flux = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum_flux -= m_flux(k,j); // total net flux
        }
        double sum_ion = 0.0;
        for (size_t k : m_kCharge) {
            sum_ion += Y(x,k,j);
        }
        // The portion of correction for ions is taken off
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += Y(x,k,j) / (1-sum_ion) * sum_flux;
        }

        // ambipolar diffusion
        double E_ambi = E(x,j);
        for (size_t k : m_kCharge) {
            double Yav = 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            double drift = rho * Yav * E_ambi
                           * m_speciesCharge[k] * m_mobility[k+m_nsp*j];
            m_flux(k,j) += drift;
        }
    }
}

void IonFlow::setSolvingStage(const size_t stage)
{
    if (stage == 1 || stage == 2 || stage == 3) {
        m_stage = stage;
    } else {
        throw CanteraError("IonFlow::updateDiffFluxes",
                    "solution phase must be set to:"
                    "1: frozenIonMethod"
                    "2: chargeNeutralityModel"
                    "3: poissonEqnMethod");
    }
}

void IonFlow::setElectricPotential(const double v1, const double v2)
{
    // This method can be used when you want to add external voltage
    m_inletVoltage = v1;
    m_outletVoltage = v2;
}

void IonFlow::setTransverseElecField(double elec_field, double elec_freq)
{
    m_elec_field = elec_field;
    m_elec_frequency = elec_freq;
}

void IonFlow::setPlasmaSourceMultiplier(double multiplier)
{
    m_plasma_multiplier = multiplier;
}

void IonFlow::setElectronTransportMultiplier(double multiplier)
{
    m_electron_multiplier = multiplier;
}

void IonFlow::enableElecHeat(bool withElecHeat)
{
    m_do_elec_heat = withElecHeat;
}

double IonFlow::getElecMobility(size_t j)
{
    return m_mobility[m_kElectron+m_nsp*j];
}

double IonFlow::getElecDiffCoeff(size_t j)
{
    return m_diff[m_kElectron+m_nsp*j];
}

double IonFlow::getElecTemperature(size_t j)
{
    return m_electronTemperature[j];
}

double IonFlow::getElecCollisionHeat(size_t j)
{
    return m_electronPower[j];
}

double IonFlow::getElecField(size_t j)
{
    return m_Eambi[j];
}

void IonFlow::setPlasmaLocation(double z1, double z2)
{
    for (size_t j = 0; j < m_points; j++) {
        m_plasmaLocation = z1;
        m_plasmaRange = z2;
        if (z(j) >= (z1-0.5*z2) && z(j) <= (z1+0.5*z2)) {
            solvePlasma(j);
        }
    }
}

void IonFlow::solvePlasma(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_plasma[i]) {
                changed = true;
            }
            m_do_plasma[i] = true;
        }
    } else {
        if (!m_do_plasma[j]) {
            changed = true;
        }
        m_do_plasma[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::solvePoissonEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = true;
        }
    } else {
        if (!m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixElectricPotential(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = false;
        }
    } else {
        if (m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_P, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::solveVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = true;
        }
    } else {
        if (!m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = false;
        }
    } else {
        if (m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::_finalize(const double* x)
{
    FreeFlame::_finalize(x);

    bool p = m_do_poisson[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!p) {
            m_fixedElecPoten[j] = phi(x, j);
        }
    }
    if (p) {
        solvePoissonEqn();
    }
    // save the velocity profile if the velocity is disabled
    bool v = m_do_velocity[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!v) {
            m_fixedVelocity[j] = u(x,j);
        }
    }
    if (v) {
        solveVelocity();
    }

    //
    bool plasma = m_do_plasma[0];
    if (plasma) {
        solvePlasma();
    }

    // for (size_t j = 0; j < m_points; j++) {
    //     if (m_do_plasma[j]) {
    //         updatePlasmaProperties(x,j);
    //     }
    // }
}

}
