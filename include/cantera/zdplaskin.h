//! @file zdplaskin.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ZDPLASKIN_H
#define CT_ZDPLASKIN_H

extern "C"
{
    void zdplaskinInit();
    void zdplaskinSetConfig(const double* atol, const double* rtol);
    void zdplaskinSetDensity(const char* cstring, const double* DENS);

    void zdplaskinSetConditions(const double* gas_temperature,
                                const double* elec_field,
                                const double* elec_frequency,
                                const double* num_density);
    double zdplaskinGetElecTemp();
    double zdplaskinGetElecDiffCoeff();
    double zdplaskinGetElecMobility(const double*);
    double zdplaskinGetElecPower(const double*);
    void zdplaskinGetPlasmaSource(double** array);
    size_t zdplaskinNSpecies();
    void zdplaskinGetSpeciesName(char* cstring, size_t* index);
}

#endif

