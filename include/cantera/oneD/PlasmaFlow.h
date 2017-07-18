//! @file PlasmaFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAFLOW_H
#define CT_PLASMAFLOW_H

#include "cantera/oneD/IonFlow.h"

using namespace std;

extern "C"
{
    void __zdplaskin_MOD_zdplaskin_init();
    void __zdplaskin_MOD_zdplaskin_set_density(const char* string,
                                               const double* DENS,
                                               const int* LDENS_CONST);
    void __zdplaskin_MOD_zdplaskin_set_conditions(const double* GAS_TEMPERATURE,
                                                  const double* REDUCED_FREQUENCY,
                                                  const double* REDUCED_FIELD,
                                                  const double* ELEC_TEMPERATURE,
                                                  const int* GAS_HEATING,
                                                  const double* SPEC_HEAT_RATIO,
                                                  const double* HEAT_SOURCE,
                                                  const int* SOFT_RESET);
    void __zdplaskin_MOD_zdplaskin_get_conditions(double* GAS_TEMPERATURE,
                                                  double* REDUCED_FREQUENCY,
                                                  double* REDUCED_FIELD,
                                                  double* ELEC_TEMPERATURE,
                                                  double* ELEC_DRIFT_VELOCITY,
                                                  double* ELEC_DIFF_COEFF,
                                                  double* ELEC_MOBILITY_N,
                                                  double* ELEC_MU_EPS_N,
                                                  double* ELEC_DIFF_EPS_N,
                                                  double* ELEC_FREQUENCY_N,
                                                  double* ELEC_POWER_N,
                                                  double* ELEC_POWER_ELASTIC_N,
                                                  double* ELEC_POWER_INELASTIC_N,
                                                  double* ELEC_EEDF);
    void __zdplaskin_MOD_zdplaskin_get_rates(double* SOURCE_TERMS,
                                             double* REACTION_RATES,
                                             double* SOURCE_TERMS_MATRIX,
                                             double* MEAN_DENSITY,
                                             double* MEAN_SOURCE_TERMS,
                                             double* MEAN_REACTION_RATES,
                                             double* MEAN_SOURCE_TERMS_MATRIX);
}
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
    void evalEEDF(size_t j, double* x, double* rsd, integer* diag, double rdt);

protected:
    virtual void updateTransport(double* x, size_t j0, size_t j1);
    vector<size_t> m_collisionSpeciesIndex;

    // ZDPlasKin wrapper
    void zdplaskin_init() {
        __zdplaskin_MOD_zdplaskin_init();
    };

    void zdplaskin_set_density(const char* speciesName,
                               const double* density,
                               const int* keepConst) {
        __zdplaskin_MOD_zdplaskin_set_density(speciesName,
                                              density,
                                              keepConst);
    };

    void zdplaskin_set_density(string speciesName, double density) {
        const double* DENS = &density;
        __zdplaskin_MOD_zdplaskin_set_density(speciesName.c_str(), DENS, NULL);
    };

    void zdplaskin_set_conditions(double gas_temperature, double reduced_field) {
        const double* GAS_TEMPERATURE = &gas_temperature;
        const double* REDUCED_FIELD = &reduced_field;
        __zdplaskin_MOD_zdplaskin_set_conditions(GAS_TEMPERATURE,
                                                 NULL, REDUCED_FIELD,
                                                 NULL, NULL, NULL, NULL, NULL);
    };

    double getElectronTemperature() {
        double* ELEC_TEMPERATURE;
        __zdplaskin_MOD_zdplaskin_get_conditions(NULL, NULL, NULL,
                                                 ELEC_TEMPERATURE,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL);
        return *ELEC_TEMPERATURE;
    };

    double getElectronDiffusionCoeff() {
        double* ELEC_DIFF_COEFF;
        __zdplaskin_MOD_zdplaskin_get_conditions(NULL, NULL, NULL, NULL,
                                                 ELEC_DIFF_COEFF,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL);
        return *ELEC_DIFF_COEFF * 1e-4;
    };

    double getElectronMobility(double number_density) {
        double* ELEC_MOBILITY_N;
        __zdplaskin_MOD_zdplaskin_get_conditions(NULL, NULL, NULL, NULL, NULL,
                                                 ELEC_MOBILITY_N,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL);
        return *ELEC_MOBILITY_N * 1e2 / number_density;
    };

    double getElectronPowerElastic(double number_density) {
        double* ELEC_POWER_ELASTIC_N;
        __zdplaskin_MOD_zdplaskin_get_conditions(NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL,
                                                 ELEC_POWER_ELASTIC_N,
                                                 NULL, NULL);
        return *ELEC_POWER_ELASTIC_N * 1e-6 * ElectronCharge * number_density;
    };

    double getElectronPowerInelastic(double number_density) {
        double* ELEC_POWER_INELASTIC_N;
        __zdplaskin_MOD_zdplaskin_get_conditions(NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 ELEC_POWER_INELASTIC_N,
                                                 NULL);
        return *ELEC_POWER_INELASTIC_N * 1e-6 * ElectronCharge * number_density;
    };

    double* getPlasmaSourceRates() {
        double* SOURCE_TERMS;
        __zdplaskin_MOD_zdplaskin_get_rates(SOURCE_TERMS,
                                            NULL, NULL, NULL, NULL,
                                            NULL, NULL);
        return SOURCE_TERMS;
    };
};

}

#endif
