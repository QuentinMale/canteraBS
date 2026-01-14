//! @file DetailedVVVTRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/DetailedVVVTRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

bool DetailedVibData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    return changed;
}

// void DetailedVibData::update(double T)
// {
//     ReactionData::update(T);
// }


DetailedVVVTRate::DetailedVVVTRate()
{
    m_B_str = "B"; //!< The string for the B parameter in VV-VT reactions
    m_C_str = "C"; //!< The string for the C parameter in VV-VT reactions
    m_D_str = "D"; //!< The string for the D parameter in VV-VT reactions
    m_scaling_str = "scaling"; //!< The string for the B parameter in VV-VT reactions
}

DetailedVVVTRate::DetailedVVVTRate(double A, double B, double C, double D, double b, double scaling)
    : ArrheniusBase(A, b, 0)
{
    writelog("On passe dans DetailedVVVTRate::DetailedVVVTRate(double A, double B, double C, double D, double b, double scaling)\n");
    m_B = B;
    m_C = C;
    m_D = D;
    m_scaling = scaling;
    
}

DetailedVVVTRate::DetailedVVVTRate(const AnyMap& node, const UnitStack& rate_units)
    : DetailedVVVTRate()
{
    // CNB: setParameters method must be implemented in TwoTempPlasmaRate
    AnyMap rate_map = node["rate-constant"].as<AnyMap>();
    UnitSystem units = node.units();

    if (rate_map.hasKey(m_B_str))
    {
        m_B = rate_map[m_B_str].asDouble();
        
    }

    if (rate_map.hasKey(m_C_str))
    {
        m_C = rate_map[m_C_str].asDouble();
        
    }

    if (rate_map.hasKey(m_D_str))
    {
        m_D = rate_map[m_D_str].asDouble();
        
    }

    if (rate_map.hasKey(m_scaling_str))
    {
        m_scaling = rate_map[m_scaling_str].asDouble();
        
    }

    setParameters(node, rate_units);
}

double DetailedVVVTRate::ddTScaledFromStruct(const DetailedVibData &shared_data) const
    {
        warn_user("DetailedVVVTRate::ddTScaledFromStruct",
                  "Temperature derivative does not consider changes of electron temperature.");
        return (m_Ea_R - m_E4_R) * shared_data.recipT * shared_data.recipT;
}

void DetailedVVVTRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // TwoTempPlasmaReaction is for a non-equilibrium plasma, and the reverse rate
    // cannot be calculated from the conventional thermochemistry.
    // @todo implement the reversible rate for non-equilibrium plasma
    if (rxn.reversible) {
        throw InputFileError("DetailedVVVTRate::setContext", rxn.input,
            "DetailedVVVTRate does not support reversible reactions");
    }
}

}
