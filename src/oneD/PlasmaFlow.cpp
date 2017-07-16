//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/PlasmaFlow.h"

using namespace std;

extern "C"
{
    double __math_MOD_getradius(double *degree);
}

namespace Cantera
{

PlasmaFlow::PlasmaFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    IonFlow(ph, nsp, points)
{
    double degree = 90.0;
    double result = __math_MOD_getradius(&degree);
    cout << result << endl;
}

void PlasmaFlow::resize(size_t components, size_t points){
    IonFlow::resize(components, points);
}

void PlasmaFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x,j0,j1);
}

void PlasmaFlow::eval(size_t jg, double* xg,
                  double* rg, integer* diagg, double rdt)
{
    IonFlow::eval(jg, xg, rg, diagg, rdt);
}

void PlasmaFlow::evalEEDF(size_t j, double* x, double* rsd, integer* diag, double rdt)
{       

}

}
