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

bool ElectronCollisionPlasmaData::update(const ThermoPhase& phase,
                                         const Kinetics& kin)
{
    const PlasmaPhase& pp = dynamic_cast<const PlasmaPhase&>(phase);

    if (pp.distributionNumber() == m_dist_number) {
        return false;
    } else {
        // update m_dist_number
        m_dist_number = pp.distributionNumber();
    }

    if (pp.levelNumber() != m_level_number) {
        levelChanged = true;
        // update m_level_number
        m_level_number = pp.levelNumber();
    } else {
        levelChanged = false;
    }

    vector<double> levels(pp.nElectronEnergyLevels());
    pp.getElectronEnergyLevels(levels.data());

    vector<double> dist(pp.nElectronEnergyLevels());
    pp.getElectronEnergyDistribution(dist.data());

    energyLevels = std::move(levels);
    distribution = std::move(dist);

    return true;
}

void ElectronCollisionPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    if (!node.hasKey("energy-levels") && !node.hasKey("cross-sections")) {
        return;
    }
    if (node.hasKey("energy-levels")) {
        m_energyLevels = node["energy-levels"].asVector<double>();
    }
    if (node.hasKey("cross-sections")) {
        m_crossSections = node["cross-sections"].asVector<double>();
    }
    if (m_energyLevels.size() != m_crossSections.size()) {
        throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
            "Energy levels and cross section must have the same length.");
    }
}

void ElectronCollisionPlasmaRate::getParameters(AnyMap& node) const {
    node["type"] = type();
    node["energy-levels"] = m_energyLevels;
    node["cross-sections"] = m_crossSections;
}

double ElectronCollisionPlasmaRate::evalFromStruct(const ElectronCollisionPlasmaData& shared_data)
{
    // Interpolate cross-sections data to the energy levels of
    // the electron energy distribution function
    if (shared_data.levelChanged) {
        m_crossSectionsInterpolated.resize(0);
        for (auto i = 0; i < shared_data.energyLevels.size(); i ++ ) {
            m_crossSectionsInterpolated.push_back(linearInterp(shared_data.energyLevels[i],
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

    // unit in kmol/m3/s
    return 0.5 * pow(2.0 * ElectronCharge / ElectronMass, 0.5) * Avogadro *
           simpson(distribution.cwiseProduct(cs_array), eps.pow(2.0));
}

void ElectronCollisionPlasmaRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    std::string electronName = getElectronSpeciesName(kin);
    if (rxn.reversible) {
        checkSpeciesCount(rxn.reactantString(), kin, electronName, rxn.input);
        checkSpeciesCount(rxn.productString(), kin, electronName, rxn.input);
        // turn on super-elastic rate
        m_enableSuperElastic = true;
    } else {
        checkSpeciesCount(rxn.reactantString(), kin, electronName, rxn.input);
    }
}

void ElectronCollisionPlasmaRate::checkSpeciesCount(const std::string& speciesString,
                                                    const Kinetics& kin,
                                                    const std::string& electronName,
                                                    AnyMap input)
{
    std::istringstream iss(speciesString);
    std::string token;
    int nElectron = 0;
    int nSpecies = 0;

    while (iss >> token) {
        if (isdigit(token[0])) {
            throw InputFileError("ElectronCollisionPlasmaRate::checkSpeciesCount", input,
                "ElectronCollisionPlasmaRate requires one electron and one molecule as species");
        }
        if (token == electronName) {
            nElectron++;
        }
        if (token != "+") {
            nSpecies++;
        }
    }

    if (nSpecies != 2) {
        throw InputFileError("ElectronCollisionPlasmaRate::checkSpeciesCount", input,
            "ElectronCollisionPlasmaRate requires exactly two species");
    }

    if (nElectron != 1) {
        throw InputFileError("ElectronCollisionPlasmaRate::checkSpeciesCount", input,
            "ElectronCollisionPlasmaRate requires one and only one electron");
    }
}

std::string ElectronCollisionPlasmaRate::getElectronSpeciesName(const Kinetics& kin) {
    if (kin.thermo().type() == "plasma") {
        return dynamic_cast<const PlasmaPhase&>(kin.thermo()).electronSpeciesName();
    } else {
        throw CanteraError("ElectronCollisionPlasmaRate::getElectronSpeciesName",
                           "ElectronCollisionPlasmaRate requires plasma phase");
    }
}

}
