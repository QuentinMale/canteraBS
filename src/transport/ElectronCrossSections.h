/**
 * @file ElectronCrossSections.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ELETRON_CROSS_SECTIONS_H
#define CT_ELETRON_CROSS_SECTIONS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Interpolate electron cross sections
/*!
 * This class provides functions that interpolate electron cross sections
 * from lxcat
 *
 * @ingroup tranprops
 */
class ElectronCrossSections
{
public:
    ElectronCrossSections() {}

private:
    doublereal fitDelta(int table, int ntstar, int degree, doublereal* c);
    static std::vector<vector_fp> electron_energy;
    static std::vector<vector_fp> effective_cross_section;
};

}

#endif
