//! @file PlasmaPressureReactor.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/PlasmaPressureReactor.h"
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

void PlasmaPressureReactor::setThermo(ThermoPhase& thermo)
{
    if (thermo.type() != "plasma") {
        throw CanteraError("PlasmaPressureReactor::setThermo",
                           "Incompatible phase type provided");
    }
    Reactor::setThermo(thermo);
    m_plasma = &dynamic_cast<PlasmaPhase&>(thermo);
    compute_disVPower();
}

void PlasmaPressureReactor::getState(double* y)
{
    if (m_plasma == 0) {
        throw CanteraError("IdealGasConstPressureReactor::getState",
                           "Error: reactor is empty.");
    }
    m_plasma->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_plasma->density() * m_vol;
    y[0] = m_mass;

    // Set the second component to the temperature
    y[1] = m_plasma->temperature();

    // set components y+2 ... y+K+1 to the vibrational energy of each species
    m_plasma->getVibrationalEnergies(y+2);

    // set components y+2+m_nspevib ... y+K+1+m_nspevib to the mass fractions of each species
    m_plasma->getMassFractions(y+2+m_nspevib);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 2 + m_nspevib);
}

// void PlasmaPressureReactor::initialize(double t0)
// {
//     IdealGasConstPressureReactor::initialize(t0);
//     printf("m_nv after base init = {}\n", m_nv);
    
//     // Number of equation in the reactor
//     // Equation for vibrational energy density is taken into account here.
//     m_nspevib = m_plasma->nsp_evib();
//     m_nv += m_nspevib;
//     disVibVPower.resize(m_nspevib);
//     RvtVPower.resize(m_nspevib);
//     recoverVibSpecies(); // get all the vibrationnal species to use this in the other functions
//     initializeStarikovskiyReading(); // initialize the reading of the yaml file for the starikovskiy model. This avoids having a heavy I/O operation at each time step.
//     printf("m_nsp = {}, m_nspevib = {}\n", m_nsp, m_nspevib);
    
// }

void PlasmaPressureReactor::initialize(double t0)
{
    IdealGasConstPressureReactor::initialize(t0);

    std::fprintf(stderr, "[CP init] base m_nv=%zu m_nsp=%zu\n", m_nv, m_nsp);
    std::fflush(stderr);

    m_nspevib = m_plasma->nsp_evib();
    std::fprintf(stderr, "[CP init] nspevib=%zu (before add)\n", m_nspevib);
    std::fflush(stderr);

    m_nv += m_nspevib;

    std::fprintf(stderr, "[CP init] after add m_nv=%zu\n", m_nv);
    std::fflush(stderr);

    disVibVPower.resize(m_nspevib);
    RvtVPower.resize(m_nspevib);
    recoverVibSpecies();
    initializeStarikovskiyReading();
}

// void PlasmaPressureReactor::updateState(double* y)
// {
//     // The components of y are [0] the total mass, [1] the total volume,
//     // [2] the temperature, [3...K+3] are the species vibrational energies,
//     // [3+m_nspevib...K+3+m_nspevib] are the mass fractions of each species,
//     // and [K+3+m_nspevib...] are the coverages of surface species on each wall.
//     m_mass = y[0];
//     m_plasma->setVibrationalEnergies(y+2);
//     m_plasma->setMassFractions_NoNorm(y+2+m_nspevib);
//     m_plasma->setState_TD(y[1], m_mass / m_vol);
//     updateConnected(true);
//     updateSurfaceState(y + m_nsp + 2 + m_nspevib);
// }

void PlasmaPressureReactor::updateState(double* y)
{
    m_mass = y[0];
    double T = y[1];

    m_plasma->setVibrationalEnergies(y + 2);
    m_plasma->setMassFractions_NoNorm(y + 2 + m_nspevib);

    // Comme IdealGasConstPressureReactor :
    m_plasma->setState_TP(T, m_pressure);
    m_vol = m_mass / m_plasma->density();

    updateConnected(false);
    updateSurfaceState(y + 2 + m_nspevib + m_nsp);
}

void PlasmaPressureReactor::eval(double time, double* LHS, double* RHS)
{
    // writelog(" **************************** Entering eval function ******************************************\n");
    double& dmdt = RHS[0]; // dm/dt (gas phase)
    double& mcpdTdt = RHS[1]; // m * c_p * dT/dt
    double* devibdt = RHS + 2; // devib/dt
    double* mdYdt = RHS + 2 + m_nspevib; // mass * dY/dt

    dmdt = 0.0;
    mcpdTdt = 0.0;


    evalWalls(time);
    m_plasma->restoreState(m_state);
    m_plasma->getPartialMolarEnthalpies(&m_hk[0]);
    const vector<double>& mw = m_plasma->molecularWeights();
    const double* Y = m_plasma->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    evalSurfaces(LHS + m_nsp + m_nspevib + 2, RHS + m_nsp + m_nspevib + 2, m_sdot.data());
    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt += mdot_surf;

    // external heat transfer
    mcpdTdt += m_Qdot;

    // gas heating from the discharge
    compute_disVPower();
    // ////////////////////////////////////////////////////////////////// DEBUG /////////////////////////////////////////////////////////////////
    // if(m_disVPower = 0){
    //     for (size_t n = 0; n < m_nsp; n++) {
    //         m_kin->
    // }

    // ////////////////////////////////////////////////////////////////// DEBUG /////////////////////////////////////////////////////////////////

    // printf(" **************************** disVPower successfulyy computed with value %f ******************************************\n", m_disVPower);
    double tot_vib_power = 0;
    // printf("nb of vib species considered : %ld \n", m_nspevib);
    if (m_nspevib > 0) {
        compute_disVibVPower();
        for (size_t n = 0; n < m_nspevib; n++){
            tot_vib_power += disVibVPower[n];
        }
    }
    // printf(" **************************** total vibrational power successfully computed with value %f ******************************************\n", tot_vib_power);
    mcpdTdt += (m_disVPower - tot_vib_power) * m_vol; // RAW FAST GAS HEATING POWER (BEFORE APPLYING PLASMA CHEMICAL SOURCE TERMS)

    // gas heating from vibrational–translational relaxation
    double tot_relax_power = 0;
    if (m_nspevib > 0) {
        compute_RvtVPower();
        
        for (size_t n = 0; n < m_nspevib; n++){
            tot_relax_power += RvtVPower[n];
        }
    }
    // printf(" **************************** total relaxation power successfuly computed with value %f ******************************************\n", tot_relax_power);
    mcpdTdt += tot_relax_power * m_vol; // SLOW GAS HEATING POWER
    

    // printf(" Entering species chemical for loop\n");
    for (size_t n = 0; n < m_nsp; n++) {
        
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        mdYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n];
        // dilution by net surface mass flux
        mdYdt[n] -= Y[n] * mdot_surf;
        //Assign left-hand side of dYdt ODE as total mass
        LHS[n+2+m_nspevib] = m_mass;
    }
    //printf("MARKER 1\n");
    for (size_t n=0; n < m_nspevib; n++){
        devibdt[n] = disVibVPower[n] - RvtVPower[n];
    }
    
    //printf("MARKER 2\n");
    // Assign left-hand side of devibdt as one
    for (size_t n = 0; n < m_nspevib; n++){
        LHS[2+n] = 1;
        //LHS[3+n] = devibdt[n];
    }
    //printf("MARKER 3\n");
    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += mdot_spec - mdot * Y[n];

            // In combination with h_in*mdot_in, flow work plus thermal
            // energy carried with the species
            mcpdTdt -= m_hk[n] / mw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        LHS[1] = m_mass * m_plasma->cp_mass();
    } else {
        RHS[1] = 0.0;
    }
    // writelog(" **************************** Exiting eval function ******************************************\n");
}

size_t PlasmaPressureReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 2 + m_nspevib;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "temperature") {
        return 1;
    } else if (nm == "evib") {
        return 2;
    } else {
        return npos;
    }
}

string PlasmaPressureReactor::componentName(size_t k) {
    if (k == 1) {
        return "temperature";
    } else if (k == 0) {
        return "mass";
    } else if (k >= 2 && k < 2 + m_nspevib) {
        return "evib";
    } else {
        return IdealGasConstPressureReactor::componentName(k - m_nspevib);
    }
}

void PlasmaPressureReactor::compute_disVPower() {
    if (m_plasma == nullptr) {
        throw CanteraError("compute_disVPower", "m_plasma is null");
    }
    if (m_plasma->E() < 1e-21){
        // If the electric field is too low, we assume no discharge power.
        m_disVPower = 0;
    }
    else{
        m_disVPower = ElectronCharge * m_plasma->nElectron()
            * m_plasma->electronMobility()
            * pow(m_plasma->E(), 2);
    }
    
    }

void PlasmaPressureReactor::compute_disVibVPower() { 
    if (m_plasma == nullptr) {
        throw CanteraError("compute_disVibVPower", "m_plasma is null");
    }
    m_kin->getNetRatesOfProgress(&m_kr[0]); // "kr"
    size_t n_vib_species = m_nspevib;
    
    for (size_t k = 0; k<n_vib_species; k++){
        disVibVPower[k] = 0;
        string vib_spec_here = vib_spec[k];
        
        for (size_t n = 0; n < m_kin->nReactions(); n++) {
            string reac_target_spec = m_plasma->getTarget(n);
            if (reac_target_spec == vib_spec_here) {
                double DUVibValue = m_plasma->getDuvib(n)*ElectronCharge; // convert to Joules
                disVibVPower[k] += DUVibValue * m_kr[n] * Avogadro; //multiply by the avogadro number to actually get a power
            }
             
        }
    }
    
     
}

std::vector<double> PlasmaPressureReactor::get_disVibVPower() {
    compute_disVibVPower();
    return disVibVPower;
}

std::vector<double> PlasmaPressureReactor::get_RvtVPower() {
    compute_RvtVPower();
    return RvtVPower;
}

std::vector<double> PlasmaPressureReactor::get_eVib() {
    if (m_plasma == nullptr) {
        throw CanteraError("get_eVib", "m_plasma is null");
    }
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


void PlasmaPressureReactor::compute_RvtVPower() {
    if (m_plasma == nullptr) {
        throw CanteraError("compute_RvtVPower", "m_plasma is null");
    }
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

double PlasmaPressureReactor::compute_TauRelax(size_t n){
    
    double tau = 0;
    string spec_name = vib_spec[n];
    // printf("Computing relaxation time for species %s\n", spec_name.c_str());
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
        throw CanteraError("PlasmaPressureReactor::compute_TauRelax",
                           "Error: species vibrational relaxation type not implemented. Please correct the YAML file or implement this species correlation.");
    }
    
    return tau;
}

double PlasmaPressureReactor::tau_millikan_white(string spec_name){ // TO BE CORRECTED OR WELL IMPLEMENTED.
    double tau = 0;
    double T = m_plasma->temperature();
    double P = m_plasma->pressure();
    double reduced_mass = 0;
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
        throw CanteraError("PlasmaPressureReactor::compute_TauRelax",
                           "Error: species vibrational relaxation time not implemented. Please correct the YAML file or implement this species correlation.");
    }
    
    epsilon_J = epsilon * ElectronCharge; // J
    theta = epsilon_J / Boltzmann; // K, the vibrational temperature of the molecule.
    
    double exponent = 5e-4 * pow(reduced_mass, 0.5) * pow(theta, 0.8) * (pow(T, -0.33) - 0.015*pow(reduced_mass, 0.25));
    double prefactor = pow(10, -8)/P;
        
    tau = prefactor * pow(10, exponent);
    
    printf("tau_millikan_%s = %e\n", spec_name.c_str(), tau);
    
    return tau;
} 

double PlasmaPressureReactor::tau_castela(string spec_name){
    double tau = 1e-11; // Almost as fast gas heating if castela is called for a species which is not N2.
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

        double p_n2 = Max(P * x_n2, 1e-16); // avoid division by zero
        double p_o2 = Max(P * x_o2, 1e-16);
        double p_o = Max(P * x_o, 1e-16);

        double tau_n2 = (exp(a_n2*(pow(T, -0.3333) - b_n2) - 18.42))*c/p_n2;
        double tau_o2 = (exp(a_o2*(pow(T, -0.3333) - b_o2) - 18.42))*c/p_o2;
        double tau_o = (exp(a_o*(pow(T, -0.3333) - b_o) - 18.42))*c/p_o;

        tau = 1/(1/tau_n2 + 1/tau_o2 + 1/tau_o);
        
    }
    // printf("tau_castela = %e\n", tau);
    return tau;
}


// Fonction pour calculer k(T) en cm3/s
double PlasmaPressureReactor::compute_k(const RelaxationEntry& entry, double T) {
    return entry.A * std::pow(T, entry.n) * std::exp(
        entry.K - entry.B / std::pow(T, 1.0 / 3.0)
                 + entry.C / std::pow(T, entry.m)
                 + entry.D / std::pow(T, entry.z)
    );
}


void PlasmaPressureReactor::readStarikovskiyRelaxYamlFile(string filename){
    // On retrouve le chemin complet à partir du nom de fichier
    std::string full_path;
    try {
        full_path = findInputFile(filename);  // cherche dans tous les chemins Cantera
    } catch (CanteraError& err) {
        throw CanteraError("PlasmaPressureReactor::readStarikovskiyRelaxYamlFile",
            "Could not find the YAML file for Starikovskiy relaxation: {}\n"
            "File requested: {}\n", err.what(), filename);
    }

    YAML::Node root = YAML::LoadFile(full_path);

    for (size_t n=0; n<m_nspevib; n++){
        string spec_name = vib_spec[n];
        // printf("Reading STARIKOVSKIY DATA YAML file for species %s\n", spec_name.c_str()); 
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
        m_data_starikovskiy.push_back(reactions);
    }
}


double PlasmaPressureReactor::tau_starikovskiy(size_t n){

    if (!starikovskiy_read) {
        if (starikovskiy_yaml_path == "init"){
            starikovskiy_yaml_path = "plasma_relax/starikovskiy_default.yaml"; // default path
            // printf("No yaml file provided for the Starikovskiy relaxation model but this model is used\n. Using default path %s\n", starikovskiy_yaml_path.c_str());
        }
        readStarikovskiyRelaxYamlFile(starikovskiy_yaml_path);
        starikovskiy_read = true;
    }

    double one_over_tau = 0;
    double T = m_plasma->temperature();
    double avogadro_per_mol = Avogadro/1000;

    // std::cout << "Reactions rates for target " << vib_spec[n] << " at T = " << T << " K:\n";
    for (const auto& r : m_data_starikovskiy[n]) {
        double k = 1e-6*compute_k(r, T); // convert to m3/s bc the result from compute_k is in cm3/s
        // std::cout << "  " << r.name << ": k = " << k << "\n";
        double x_partner = m_plasma->moleFraction(r.name); // loop over all the species and moleFraction returns 0 if the species is not in the phase
        double mixture_molar_density = 1000*m_plasma->molarDensity(); // in cantera, the density in kmol/m^3 so we need to convert to mol/m^3 to get things right
        one_over_tau += k * x_partner * mixture_molar_density * avogadro_per_mol;
    }
    double tau = 1/one_over_tau;
    // std::cout << "Starikovskiy relaxation time for " << vib_spec[n] << ": " << tau << "\n";

    return tau;
} 

void  PlasmaPressureReactor::recoverVibSpecies(){
    vib_spec = m_plasma->getVibSpecies();
    // printf("Vibrational species recovery:\n");
    if (vib_spec.size() == 0){
        // printf("    No vibrational species found\n");
    }
    for (size_t n = 0; n < vib_spec.size(); n++){
        // printf("Vibrational species %ld: %s\n", n, vib_spec[n].c_str());
    }
}

void PlasmaPressureReactor::setVibRelaxType(string relax_type_name){
    relax_type = relax_type_name;
    // printf("Relaxation type set to %s\n", relax_type.c_str());
    }

string PlasmaPressureReactor::getVibRelaxType(){
    return relax_type;
}

double PlasmaPressureReactor::getVibConstantModelTauRelax(){
    return tau_relax_constant_model;
}

void PlasmaPressureReactor::setVibConstantModelTauRelax(double tau_to_set){
    tau_relax_constant_model = tau_to_set;
    // printf("Relaxation time constant model set to %f\n", tau_relax_constant_model);
    }

double PlasmaPressureReactor::Max(double a, double b){
    if (a>b) {
        return a;
    } else{
        return b;
    }
    }

}




