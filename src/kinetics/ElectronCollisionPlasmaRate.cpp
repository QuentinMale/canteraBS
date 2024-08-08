//! @file ElectronCollisionPlasmaRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/numerics/funcs.h"

namespace Cantera
{

ElectronCollisionPlasmaData::ElectronCollisionPlasmaData()
{
    energyLevels.assign(1, 0.0);
    distribution.assign(1, 0.0);
}

bool ElectronCollisionPlasmaData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    const PlasmaPhase& pp = dynamic_cast<const PlasmaPhase&>(phase);

    // The distribution number dictates whether the rate should be updated.
    // Three scenarios involving changes of the distribution number:
    // 1. Change of the electron energy levels
    // 2. Change of the electron energy distribution
    // 3. Combined changes of one and two
    if (pp.distributionNumber() == m_dist_number) {
        return false;
    }

    // Update distribution
    m_dist_number = pp.distributionNumber();
    distribution.resize(pp.nElectronEnergyLevels());
    pp.getElectronEnergyDistribution(distribution.data());

    // Update energy levels
    // levelChanged = pp.levelNumber() != m_level_number;
    // if (levelChanged) {
        m_level_number = pp.levelNumber();
        energyLevels.resize(pp.nElectronEnergyLevels());
        pp.getElectronEnergyLevels(energyLevels.data());
    // }

    return true;
}

void ElectronCollisionPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);

    // Two options are available to define the cross sections of the Electron Collision Plasma Rate:
    // 1) The data are retrieved from a set of cross-section provided in the yaml file, in that case
    //    the kind of the process, as well as target and product of the process must be provided in
    //    the reaction and must match with one of the loaded cross-section. In that case, the data are 
    //    loaded in PlasmaPhase::initThermo. Recommended when the EEDF is solved
    // 2) The energy levels and cross sections are directly given in the reaction. This can be usefull
    //    when the EEDF is not solved (no set of cross-sections loaded) or if you want to include cross-
    //    section that is not part of the cross-section dataset or if you want to force some cross-sections.

    // Method 1:
    // Retrieve the associated process data (must match with one loaded cross-section)
    // Mandatory if energy levels and cross-sections are note given in the reaction
    if (node.hasKey("kind"))
    {
        m_kind = node["kind"].asString();
    }
    if (node.hasKey("target"))
    {
        m_target = node["target"].asString();
    }
    if (node.hasKey("product"))
    {
        m_product = node["product"].asString();
    }

    // Method 2:
    if (node.hasKey("energy-levels") && node.hasKey("cross-sections"))
    {   
        writelog("WARNING: cross section is found in the reaction and imposed!\n");
        if (node.hasKey("energy-levels"))
        {
            m_energyLevels = node["energy-levels"].asVector<double>();
        }
        if (node.hasKey("cross-sections"))
        {
            m_crossSections = node["cross-sections"].asVector<double>();
        }
        if (m_energyLevels.size() != m_crossSections.size())
        {
            throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
                                "Energy levels and cross section must have the same length.");
        }
        cs_ok = true;
    }
}

void ElectronCollisionPlasmaRate::getParameters(AnyMap& node) const {
    node["type"] = type();
    node["energy-levels"] = m_energyLevels;
    node["cross-sections"] = m_crossSections;
}

double ElectronCollisionPlasmaRate::evalFromStruct(
    const ElectronCollisionPlasmaData& shared_data)
{
    // Interpolate cross-sections data to the energy levels of
    // the electron energy distribution function when the interpolated
    // cross section is empty
    // Note that the PlasmaPhase should handle the interpolated cross sections
    // for all ElectronCollisionPlasmaRate reactions
    if (m_crossSectionsInterpolated.size() == 0) {
        for (double level : shared_data.energyLevels) {
            m_crossSectionsInterpolated.push_back(linearInterp(level,
                                                  m_energyLevels, m_crossSections));
        }
    }

    // Map cross sections to Eigen::ArrayXd
    auto cs_array = Eigen::Map<const Eigen::ArrayXd>(
        m_crossSectionsInterpolated.data(), m_crossSectionsInterpolated.size()
    );

    // Map energyLevels in Eigen::ArrayXd
    auto eps = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.energyLevels.data(), shared_data.energyLevels.size()
    );

    // Map energyLevels in Eigen::ArrayXd
    auto distribution = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.distribution.data(), shared_data.distribution.size()
    );

    // Integrate reaction rate (unit in kmol/m3/s)
    string quadratureMethod = "simpson";
    Eigen::VectorXd y(distribution.size());
    for (size_t i = 0; i < distribution.size(); i++)
    {
        y[i] = cs_array[i] * eps[i] * distribution[i];
    }
    return pow(2.0 * ElectronCharge / ElectronMass, 0.5) * numericalQuadrature(quadratureMethod, y, eps) * Avogadro;

    
    // return 0.5 *
    //         pow(2.0 * ElectronCharge / ElectronMass, 0.5) * Avogadro *
    //         simpson(distribution.cwiseProduct(cs_array), eps.pow(2.0));
}

void ElectronCollisionPlasmaRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // ElectronCollisionPlasmaReaction is for a non-equilibrium plasma, and the reverse rate
    // cannot be calculated from the conventional thermochemistry.
    // @todo implement the reversible rate for non-equilibrium plasma
    if (rxn.reversible) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate does not support reversible reactions");
    }

    // get electron species name
    string electronName;
    if (kin.thermo().type() == "plasma") {
        electronName = dynamic_cast<const PlasmaPhase&>(kin.thermo()).electronSpeciesName();
    } else {
        throw CanteraError("ElectronCollisionPlasmaRate::setContext",
                           "ElectronCollisionPlasmaRate requires plasma phase");
    }

    // Number of reactants needs to be two
    if (rxn.reactants.size() != 2) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires exactly two reactants");
    }

    // Must have only one electron
    // @todo add electron-electron collision rate
    if (rxn.reactants.at(electronName) != 1) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires one and only one electron");
    }
}

}
