//! @file PlasmaReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAREACTOR_H
#define CT_PLASMAREACTOR_H

#include "IdealGasReactor.h"
#include "cantera/thermo/PlasmaPhase.h"

namespace Cantera
{

/**
 * Class PlasmaReactor is a class for stirred reactors that ...
 *
 *
 * @ingroup reactorGroup
 */
class PlasmaReactor : public IdealGasReactor
{
public:
    using IdealGasReactor::IdealGasReactor; // inherit constructors

    string type() const override {
        return "PlasmaReactor";
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

    // Structure pour stocker les coefficients d'une réaction
    struct RelaxationEntry {
        std::string name;
        std::string target;
        double A, n, K, B, C, m, D, z;
    };

    double compute_k(const RelaxationEntry& entry, double T);

    string stari_yaml_path;  // attribut à ajouter dans ton .h

    void setStariYamlPath(string path) {
        stari_yaml_path = path;
    }

    string getStariYamlPath() {
        return stari_yaml_path;
}

    void readStariRelaxYamlFile(string filename);

    void initializeStariReading(){
        stari_read = false;
    }


protected:
    void setThermo(ThermoPhase& thermo) override;

    double m_dis_vol; //!< Discharge volume

    double m_disVPower; //!< Volumetric discharge power

    std::vector<double> disVibVPower; //!< Volumetric discharge power going into vibrational excitation

    std::vector<double> RvtVPower; // Vibrational energy relaxation into heat

    size_t m_nspevib; //!< Number of species with vibrational excitation

    std::vector<std::string> vib_spec; 

    PlasmaPhase* m_plasma = nullptr;

    string relax_type;

    double tau_relax_constant_model;

    std::vector<std::vector<RelaxationEntry>> m_data_stari;

    bool stari_read;


};
}

#endif
