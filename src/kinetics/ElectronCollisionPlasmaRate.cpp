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
    EnergyLevels.resize(1);
    EnergyLevels << 0.0;
    Distribution.resize(1);
    Distribution << 0.0;
}

bool ElectronCollisionPlasmaData::update(const ThermoPhase& phase,
                                         const Kinetics& kin)
{
    const PlasmaPhase& pp = dynamic_cast<const PlasmaPhase&>(phase);

    Eigen::ArrayXd levels(pp.nElectronEnergyLevels());
    pp.getElectronEnergyLevels(levels.data());

    Eigen::ArrayXd dist(pp.nElectronEnergyLevels());
    pp.getElectronEnergyDistribution(dist.data());

    bool changed = false;

    if (levels.size() != EnergyLevels.size() || (levels - EnergyLevels).any() != 0.0) {
        EnergyLevels = std::move(levels);
        changed = true;
    }

    if (dist.size() != Distribution.size() || (dist - Distribution).any() != 0.0) {
        Distribution = std::move(dist);
        changed = true;
    }

    return changed;
}

ElectronCollisionPlasmaRate::ElectronCollisionPlasmaRate()
{
}

void ElectronCollisionPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    if (node.hasKey("energy-levels")) {
        vector<double> levels = node["energy-levels"].asVector<double>();
        m_energyLevels = Eigen::Map<const Eigen::ArrayXd>(levels.data(), levels.size());
    } else {
        m_energyLevels = NAN;
        return;
    }
    if (node.hasKey("cross-section")) {
        vector<double> crossSection = node["cross-section"].asVector<double>();
        m_crossSection = Eigen::Map<const Eigen::ArrayXd>(crossSection.data(), crossSection.size());
    } else {
        m_crossSection = NAN;
        return;
    }
    if (m_energyLevels.size() != m_crossSection.size()) {
        throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
            "Energy levels and cross section must have the same length.");
    }
}

void ElectronCollisionPlasmaRate::getParameters(AnyMap& node) const {
    node["type"] = type();
    vector<double> energyLevels(m_energyLevels.size());
    Eigen::Map<Eigen::ArrayXd>(energyLevels.data(), energyLevels.size()) = m_energyLevels;
    node["energy-levels"] = energyLevels;
    vector<double> crossSection(m_crossSection.size());
    Eigen::Map<Eigen::ArrayXd>(crossSection.data(), crossSection.size()) = m_crossSection;
    node["cross-section"] = crossSection;
}

double ElectronCollisionPlasmaRate::evalFromStruct(const ElectronCollisionPlasmaData& shared_data) const
{
    Eigen::ArrayXd eps2 = shared_data.EnergyLevels.pow(2.0);
    Eigen::ArrayXd cross_section(shared_data.EnergyLevels.size());
    linearInterp(shared_data.EnergyLevels, m_energyLevels, m_crossSection, cross_section);
    return 0.5 * pow(2.0 * ElectronCharge / ElectronMass, 0.5) *
           simpson(shared_data.Distribution.cwiseProduct(cross_section), eps2);
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
    if (rxn.reactants.size() != 2) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires exactly two reactants");
    }
    // get electron species name
    string electronName;
    if (kin.thermo().type() == "plasma") {
        electronName = dynamic_cast<const PlasmaPhase&>(kin.thermo()).electronSpeciesName();
    } else {
        throw CanteraError("ElectronCollisionPlasmaRate::setContext",
                           "ElectronCollisionPlasmaRate requires plasma phase");
    }
    // find electron in reactants
    if (rxn.reactants.find(electronName) == rxn.reactants.end()) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires at least an electron species as "
            "one reactant");
    }
}

}
