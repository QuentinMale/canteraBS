//! @file PlasmaPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMA_PRESSURE_REACTOR_H
#define CT_PLASMA_PRESSURE_REACTOR_H

#include "IdealGasConstPressureReactor.h"
#include "cantera/thermo/PlasmaPhase.h"

namespace Cantera
{

/**
 * Class PlasmaReactor is a class for stirred reactors that ...
 *
 *
 * @ingroup reactorGroup
 */
class PlasmaPressureReactor : public IdealGasConstPressureReactor
{
public:
    using IdealGasConstPressureReactor::IdealGasConstPressureReactor; // inherit constructors

    string type() const override {
        return "PlasmaPressureReactor";
    }

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void updateState(double* y) override;

    void eval(double t, double* LHS, double* RHS) override;

    //! Set/Get discharge volume
    void setDisVol(double dis_vol) {
        m_dis_vol = dis_vol;
    }
    double disVol() const {
        return m_dis_vol;
    }

    //! Get discharge volumetric power
    //CQM may not be up to date
    double disVPower() const{
        return m_disVPower;
    }

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

    void compute_disVPower();

    void compute_disVibVPower();

    void compute_RvtVPower();

    // double compute_TauRelax_N2();

    // double compute_TauRelax_O2();

    double  compute_TauRelax(size_t n);

    std::vector<double> get_disVibVPower();

    std::vector<double> get_RvtVPower();

    std::vector<double> get_eVib();

    void recoverVibSpecies();

    void setVibRelaxType(string relax_type_name);

    string getVibRelaxType();

    double getVibConstantModelTauRelax();

    void setVibConstantModelTauRelax(double tau_to_set);

    double tau_millikan_white(string spec_name);

    double tau_castela(string spec_name);

    double tau_starikovskiy(size_t n);

    // a structure to store the relaxation time data for the starikovski model. The data will be read from a yaml file provided by the user.
    struct RelaxationEntry {
        std::string name;
        std::string target;
        double A, n, K, B, C, m, D, z;
    };

    double compute_k(const RelaxationEntry& entry, double T);

    string starikovskiy_yaml_path = "init";  // the path to the yaml file containing the relaxation data for the Starikovskiy model provided by the user

    void setStarikovskiyYamlPath(string path) {
        starikovskiy_yaml_path = path;
    }

    string getStarikovskiyYamlPath() {
        return starikovskiy_yaml_path;
}

    void readStarikovskiyRelaxYamlFile(string filename);

    void initializeStarikovskiyReading(){
        starikovskiy_read = false;
    }

    double Max(double a, double b);



protected:
    void setThermo(ThermoPhase& thermo) override;

    double m_dis_vol = 1; //!< Discharge volume

    double m_disVPower = 0; //!< Volumetric discharge power

    std::vector<double> disVibVPower; //!< Volumetric discharge power going into vibrational excitation

    std::vector<double> RvtVPower; // Vibrational energy relaxation into heat

    size_t m_nspevib; //!< Number of species with vibrational excitation

    std::vector<std::string> vib_spec; // a vector to store the names of vibrational species

    PlasmaPhase* m_plasma = nullptr; // pointer to the plasma phase

    string relax_type = "Castela"; // relaxation type to be chosen by the user. It will be Castela, Starikovski, Constant or MillikanandWhite

    double tau_relax_constant_model = 1e-4; // relaxation time for the constant model

    std::vector<std::vector<RelaxationEntry>> m_data_starikovskiy; // relaxation data input from the relaxation yaml file provided by the user. Used by the Starikovskiy model

    bool starikovskiy_read = false; // boolean to check if the yaml file has been read or not, to avoid reading it several times

    Kinetics* m_kinetics;


};
}

#endif
