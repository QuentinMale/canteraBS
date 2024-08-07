//! @file PlasmaReactor.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/PlasmaReactor.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void PlasmaReactor::setThermo(ThermoPhase& thermo)
{
    if (thermo.type() != "plasma") {
        throw CanteraError("PlasmaReactor::setThermo",
                           "Incompatible phase type provided");
    }
    Reactor::setThermo(thermo);
}



}
