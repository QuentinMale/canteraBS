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
#include "cantera/base/ct_defs.h"  // contient findInputFile()


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <yaml-cpp/yaml.h>

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
    disVibVPower.resize(m_nspevib);
    RvtVPower.resize(m_nspevib);
    recoverVibSpecies(); // get all the vibrationnal species to use this in the other functions
    initializeStariReading(); // initialize the reading of the yaml file for the starikovskiy model. This avoids having a heavy I/O operation at each time step.
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
    //printf(" **************************** Entering eval function ******************************************\n");
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
    double tot_vib_power = 0;
    // printf("nb of vib species considered : %ld \n", m_nspevib);
    if (m_nspevib > 0) {
        compute_disVibVPower();
        for (size_t n = 0; n < m_nspevib; n++){
            tot_vib_power += disVibVPower[n];
        }
    }
    printf(" **************************** total vibrational power successfully computed with value %f ******************************************\n", tot_vib_power);
    mcvdTdt += (m_disVPower - tot_vib_power) * m_vol; // FAST GAS HEATING POWER

    // gas heating from vibrational–translational relaxation
    double tot_relax_power = 0;
    if (m_nspevib > 0) {
        compute_RvtVPower();
        
        for (size_t n = 0; n < m_nspevib; n++){
            tot_relax_power += RvtVPower[n];
        }
    }
    printf(" **************************** total relaxation power successfuly computed with value %f ******************************************\n", tot_relax_power);
    mcvdTdt += tot_relax_power * m_vol; // SLOW GAS HEATING POWER
    

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
    //printf("MARKER 1\n");
    for (size_t n=0; n < m_nspevib; n++){
        devibdt[n] = disVibVPower[n] - RvtVPower[n];
    }
    
    //printf("MARKER 2\n");
    // Assign left-hand side of devibdt as one
    for (size_t n = 0; n < m_nspevib; n++){
        LHS[3+n] = 1;
        //LHS[3+n] = devibdt[n];
    }
    //printf("MARKER 3\n");
    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

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

    RHS[1] = m_vdot;
    if (m_energy) {
        LHS[2] = m_mass * m_plasma->cv_mass();
    } else {
        RHS[2] = 0;
    }
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

void PlasmaReactor::compute_disVibVPower() { 
;
    m_kin->getNetRatesOfProgress(&m_kr[0]); // "kr"
    size_t n_vib_species = m_nspevib;
    
    for (size_t k = 0; k<n_vib_species; k++){
        disVibVPower[k] = 0;
        string vib_spec_here = vib_spec[k];
        
        for (size_t n = 0; n < m_kin->nReactions(); n++) {
            string reac_target_spec = m_plasma->getTarget(n);
            if (reac_target_spec == vib_spec_here) {
                double DUVibValue = m_plasma->getDuvib(n)*1.6e-19; // convert to Joules
                disVibVPower[k] += DUVibValue * m_kr[n] * 6.02e26; //multiply by the avogadro number to actually get a power
            }
             
        }
    }
    
     
}

std::vector<double> PlasmaReactor::get_disVibVPower() {
    compute_disVibVPower();
    return disVibVPower;
}

std::vector<double> PlasmaReactor::get_RvtVPower() {
    compute_RvtVPower();
    return RvtVPower;
}

std::vector<double> PlasmaReactor::get_eVib() {
    size_t n_vib_species = m_nspevib;
    std::vector<double> to_return(n_vib_species);

    double* evib_array = new double[n_vib_species];
    m_plasma->getVibrationalEnergies(evib_array);

    for (size_t n = 0; n < n_vib_species; ++n) {
        to_return[n] = evib_array[n];
    }

    delete[] evib_array;
    return to_return;
}


void PlasmaReactor::compute_RvtVPower() {
    size_t n_vib_species = m_nspevib;

    double* evib_array = new double[n_vib_species];
    
    m_plasma->getVibrationalEnergies(evib_array);
    
    for (size_t n=0; n<n_vib_species; n++){
        RvtVPower[n] = 0;
        double tau = compute_TauRelax(n);
        RvtVPower[n] = evib_array[n]/tau;
    }

    delete[] evib_array;
}

double PlasmaReactor::compute_TauRelax(size_t n){
    
    double tau = 0;
    string spec_name = vib_spec[n];
    printf("Computing relaxation time for species %s\n", spec_name.c_str());
    if (relax_type == "Millikan&White"){
        tau = tau_millikan_white(spec_name);
    }
    else if (relax_type == "Castela"){
        tau = tau_castela(spec_name);
    }
    else if (relax_type == "Constant"){
        tau = tau_relax_constant_model;
    }
    else if (relax_type == "Starikovskiy"){
        tau = tau_starikovskiy(n);
    }
    
    else{
        throw CanteraError("PlasmaReactor::compute_TauRelax",
                           "Error: species vibrational relaxation type not implemented. Please correct the YAML file or implement this species correlation.");
    }
    
    return tau;
}

double PlasmaReactor::tau_millikan_white(string spec_name){
    double tau = 0;
    double T = m_plasma->temperature();
    double P = m_plasma->pressure();
    double reduced_mass = 0;
    double boltzmann_cst = 1.38e-23; // J/K
    double epsilon = 0;
    double epsilon_J = 0;
    double theta = 0;
    
    if (spec_name == "N2"){
        reduced_mass = 1.16e-26; //kg
        epsilon = 1.21; // eV
    }
    else if (spec_name == "O2"){
        reduced_mass = 1.33e-26; //kg
        epsilon = 0.41; // eV
    }
    
    else{
        throw CanteraError("PlasmaReactor::compute_TauRelax",
                           "Error: species vibrational relaxation time not implemented. Please correct the YAML file or implement this species correlation.");
    }
    
    epsilon_J = epsilon * 1.6e-19; // J
    theta = epsilon_J / boltzmann_cst; // K, the vibrational temperature of the molecule.
    
    double exponent = 5e-4 * pow(reduced_mass, 0.5) * pow(theta, 0.8) * (pow(T, -0.33) - 0.015*pow(reduced_mass, 0.25));
    double prefactor = pow(10, -8)/P;
        
    tau = prefactor * pow(10, exponent);
    
    printf("tau_millikan_%s = %e\n", spec_name.c_str(), tau);
    
    return tau;
} 

double PlasmaReactor::tau_castela(string spec_name){
    double tau = 1e-12;
    if (spec_name == "N2") {
        double a_n2 = 221.0;
        double b_n2 = 0.029;
        double a_o2 = 229.0;
        double b_o2 = 0.0295;
        double a_o = 72.4; 
        double b_o = 0.015;
        
        double c = 101325; //Pa.s

        double T = m_plasma->temperature();
        double P = m_plasma->pressure();
        double x_n2 = m_plasma->moleFraction("N2");
        double x_o2 = m_plasma->moleFraction("O2");
        double x_o = m_plasma->moleFraction("O");

        double p_n2 = P * x_n2;
        double p_o2 = P * x_o2;
        double p_o = P * x_o;

        double tau_n2 = (exp(a_n2*(pow(T, -0.3333) - b_n2) - 18.42))*c/p_n2;
        double tau_o2 = (exp(a_o2*(pow(T, -0.3333) - b_o2) - 18.42))*c/p_o2;
        double tau_o = (exp(a_o*(pow(T, -0.3333) - b_o) - 18.42))*c/p_o;

        tau = 1/(1/tau_n2 + 1/tau_o2 + 1/tau_o);
        
    }
    printf("tau_castela = %e\n", tau);
    return tau;
}


// Fonction pour calculer k(T) en cm3/s
double PlasmaReactor::compute_k(const RelaxationEntry& entry, double T) {
    return entry.A * std::pow(T, entry.n) * std::exp(
        entry.K - entry.B / std::pow(T, 1.0 / 3.0)
                 + entry.C / std::pow(T, entry.m)
                 + entry.D / std::pow(T, entry.z)
    );
}


void PlasmaReactor::readStariRelaxYamlFile(string filename){
    // On retrouve le chemin complet à partir du nom de fichier
    std::string full_path;
    try {
        full_path = findInputFile(filename);  // cherche dans tous les chemins Cantera
    } catch (CanteraError& err) {
        throw CanteraError("PlasmaReactor::readStariRelaxYamlFile",
            "Could not find the YAML file for Starikovskiy relaxation: {}\n"
            "File requested: {}\n", err.what(), filename);
    }

    YAML::Node root = YAML::LoadFile(full_path);

    for (size_t n=0; n<m_nspevib; n++){
        string spec_name = vib_spec[n];
        printf("Reading STARI YAML file for species %s\n", spec_name.c_str()); 
        std::string key = spec_name + "_relaxations";

        std::vector<RelaxationEntry> reactions;
        for (const auto& node : root[key]) {
            RelaxationEntry r;
            r.name = node["name"].as<std::string>();
            r.target = node["target"].as<std::string>();
            r.A = node["A"].as<double>();
            r.n = node["n"].as<double>();
            r.K = node["K"].as<double>();
            r.B = node["B"].as<double>();
            r.C = node["C"].as<double>();
            r.m = node["m"].as<double>();
            r.D = node["D"].as<double>();
            r.z = node["z"].as<double>();
            reactions.push_back(r);
        }
        m_data_stari.push_back(reactions);
    }
}


double PlasmaReactor::tau_starikovskiy(size_t n){

    if (!stari_read) {
        if (stari_yaml_path == "init"){
            stari_yaml_path = "plasma_relax/stari_default.yaml"; // default path
            printf("No yaml file provided for the Starikovskiy relaxation model but this model is used\n. Using default path %s\n", stari_yaml_path.c_str());
        }
        readStariRelaxYamlFile(stari_yaml_path);
        stari_read = true;
    }

    double one_over_tau = 0;
    double T = m_plasma->temperature();
    double avogadro = 6.022e23; // in m3/mol

    // std::cout << "Reactions rates for target " << vib_spec[n] << " at T = " << T << " K:\n";
    for (const auto& r : m_data_stari[n]) {
        double k = 1e-6*compute_k(r, T); // convert to m3/s bc the result from compute_k is in cm3/s
        // std::cout << "  " << r.name << ": k = " << k << "\n";
        double x_partner = m_plasma->moleFraction(r.name);
        double mixture_molar_density = 1000*m_plasma->molarDensity(); // in cantera, the density in kmol/m^3 so we need to convert to mol/m^3 to get things right
        one_over_tau += k * x_partner * mixture_molar_density * avogadro;
    }
    double tau = 1/one_over_tau;
    std::cout << "Starikovskiy relaxation time for " << vib_spec[n] << ": " << tau << "\n";

    return tau;
} 

void  PlasmaReactor::recoverVibSpecies(){
    vib_spec = m_plasma->getVibSpecies();
    printf("Vibrational species recovery:\n");
    if (vib_spec.size() == 0){
        printf("    No vibrational species found\n");
    }
    for (size_t n = 0; n < vib_spec.size(); n++){
        printf("Vibrational species %ld: %s\n", n, vib_spec[n].c_str());
    }
}

void PlasmaReactor::setVibRelaxType(string relax_type_name){
    relax_type = relax_type_name;
    printf("Relaxation type set to %s\n", relax_type.c_str());}

string PlasmaReactor::getVibRelaxType(){
    return relax_type;
}

double PlasmaReactor::getVibConstantModelTauRelax(){
    return tau_relax_constant_model;
}

void PlasmaReactor::setVibConstantModelTauRelax(double tau_to_set){
    tau_relax_constant_model = tau_to_set;
    printf("Relaxation time constant model set to %f\n", tau_relax_constant_model);}

}




