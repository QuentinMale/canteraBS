//! @file PlasmaReactor.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/PlasmaReactor.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"

namespace Cantera
{

void PlasmaReactor::setThermo(ThermoPhase& thermo)
{
    if (thermo.type() != "plasma") {
        throw CanteraError("PlasmaReactor::setThermo",
                           "Incompatible phase type provided");
    }
    Reactor::setThermo(thermo);
    m_plasma = &dynamic_cast<PlasmaPhase&>(thermo);
    compute_disVPower();
}

void PlasmaReactor::getState(double* y)
{
    if (m_plasma == 0) {
        throw CanteraError("IdealGasReactor::getState",
                           "Error: reactor is empty.");
    }
    m_plasma->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_plasma->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // Set the third component to the temperature
    y[2] = m_plasma->temperature();

    // set components y+3 ... y+K+2 to the vibrational energy of each species
    m_plasma->getVibrationalEnergies(y+3);

    // set components y+3+m_nspevib ... y+K+2+m_nspevib to the mass fractions of each species
    m_plasma->getMassFractions(y+3+m_nspevib);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 3 + m_nspevib);
}

Kinetics* m_kinetics;

void PlasmaReactor::initialize(double t0)
{
    IdealGasReactor::initialize(t0);

    // Number of equation in the reactor
    // Equation for vibrational energy density is taken into account here.
    m_nspevib = m_plasma->nsp_evib();
    m_nv += m_nspevib;
}

void PlasmaReactor::updateState(double* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the species vibrational energies,
    // [3+m_nspevib...K+3+m_nspevib] are the mass fractions of each species,
    // and [K+3+m_nspevib...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
    m_plasma->setVibrationalEnergies(y+3);
    m_plasma->setMassFractions_NoNorm(y+3+m_nspevib);
    m_plasma->setState_TD(y[2], m_mass / m_vol);
    updateConnected(true);
    updateSurfaceState(y + m_nsp + 3 + m_nspevib);
}

void PlasmaReactor::eval(double time, double* LHS, double* RHS)
{
    printf(" **************************** Entering eval function ******************************************\n");
    double& dmdt = RHS[0]; // dm/dt (gas phase)
    double& mcvdTdt = RHS[2]; // m * c_v * dT/dt
    double* devibdt = RHS + 3; // devib/dt
    double* mdYdt = RHS + 3 + m_nspevib; // mass * dY/dt

    evalWalls(time);
    m_plasma->restoreState(m_state);
    m_plasma->getPartialMolarIntEnergies(&m_uk[0]);
    const vector<double>& mw = m_plasma->molecularWeights();
    const double* Y = m_plasma->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    evalSurfaces(LHS + m_nsp + m_nspevib + 3, RHS + m_nsp + m_nspevib + 3, m_sdot.data());
    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt += mdot_surf;

    // compression work and external heat transfer
    mcvdTdt += - m_pressure * m_vdot + m_Qdot;

    // gas heating from the discharge
    compute_disVPower();
    printf(" **************************** disVPower successfulyy computed with value %f ******************************************\n", m_disVPower);
    compute_disVibVPower();
    printf(" **************************** vibrationnal power successfully computed with value %f ******************************************\n", disVibVPower);
    mcvdTdt += (m_disVPower - disVibVPower) * m_vol; // FAST GAS HEATING POWER

    // gas heating from vibrational–translational relaxation
    printf(" **************************** Attempting to compute relaxation power  ******************************************\n");
    compute_RvtVPower();
    printf(" **************************** relaxation power successfuly computed with value %f ******************************************\n", RvtVPower);
    mcvdTdt += RvtVPower * m_vol; // SLOW GAS HEATING POWER

    printf(" Entering species chemical for loop\n");
    for (size_t n = 0; n < m_nsp; n++) {
        
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        mdYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n];
        // dilution by net surface mass flux
        mdYdt[n] -= Y[n] * mdot_surf;
        //Assign left-hand side of dYdt ODE as total mass
        LHS[n+3+m_nspevib] = m_mass;
    }
    printf("MARKER 1\n");
    devibdt[0] = disVibVPower - RvtVPower;
    printf("MARKER 2\n");
    // Assign left-hand side of devibdt as one
    for (size_t n = 0; n < m_nspevib; n++){
        LHS[3+n] = 1;
        //LHS[3+n] = devibdt[n];
    }
    printf("MARKER 3\n");
    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }
    printf("MARKER 4\n");

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        mcvdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += mdot_spec - mdot * Y[n];

            // In combination with h_in*mdot_in, flow work plus thermal
            // energy carried with the species
            mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
        }
    }
    printf("MARKER 5\n");

    RHS[1] = m_vdot;
    if (m_energy) {
        LHS[2] = m_mass * m_plasma->cv_mass();
    } else {
        RHS[2] = 0;
    }
    printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ EXITING eval function ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
}

size_t PlasmaReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 3 + m_nspevib;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else if (nm == "temperature") {
        return 2;
    } else if (nm == "evib") {
        return 3;
    } else {
        return npos;
    }
}

string PlasmaReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else if (k == 0) {
        return "mass";
    } else if (k == 1) {
        return "volume";
    } else if (k >= 3 && k < 3 + m_nspevib) {
        return "evib";
    } else {
        return Reactor::componentName(k-m_nspevib);
    }
}

void PlasmaReactor::compute_disVPower() {
    m_disVPower = ElectronCharge * m_plasma->nElectron()
            * m_plasma->electronMobility()
            * pow(m_plasma->E(), 2);
    }

void PlasmaReactor::compute_disVibVPower() { // TO IMPLEMENT
    disVibVPower = 0;
    m_kin->getNetRatesOfProgress(&m_kr[0]); // "kr"

    for (size_t n = 0; n < m_kin->nReactions(); n++) {
        
        double DUVibValue = m_plasma->getDuvib(n)*1.6e-19; // convert to Joules
        disVibVPower += DUVibValue * m_kr[n] * 6.02e26; //multiply by the avogadro number to actually get a power 
    }
    
}

double PlasmaReactor::get_disVibVPower() {
    compute_disVibVPower();
    return disVibVPower;
}

double PlasmaReactor::get_RvtVPower() {
    compute_RvtVPower();
    return RvtVPower;
}

double PlasmaReactor::get_eVib() {

    size_t n_vib_species = m_plasma->nsp_evib();
    double* evib_array = new double[n_vib_species];
    m_plasma->getVibrationalEnergies(evib_array);

    double eVib = evib_array[0];

    delete[] evib_array;

    return eVib;
}

void PlasmaReactor::compute_RvtVPower() { // TO IMPLEMENT
    RvtVPower = 0;

    size_t n_vib_species = m_plasma->nsp_evib();

    double* evib_array = new double[n_vib_species];
    // printf("Hey");
    m_plasma->getVibrationalEnergies(evib_array);
    // printf("Hoy");
    compute_TauRelaxN2();
     // printf("Hiyaaa");
    RvtVPower = evib_array[0]/tau_relax_vib_N2;

    delete[] evib_array;

}
    
void PlasmaReactor::compute_TauRelaxN2() {

    tau_relax_vib_N2 = 0.0001;

}

}
