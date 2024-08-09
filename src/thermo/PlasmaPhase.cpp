//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"
#include "cantera/numerics/funcs.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include <boost/polymorphic_pointer_cast.hpp>
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"

namespace Cantera {

PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
{
    initialize();

    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyLevels = Eigen::ArrayXd::LinSpaced(m_nPoints, 0.0, 1.0);

    // initial electron temperature
    setElectronTemperature(temperature());

    //CQM TODO set m_nspevib
    m_evib.resize(m_nspevib);

}

void PlasmaPhase::initialize()
{
    m_ncs = 0;
    m_f0_ok = false;
    m_EN = 0.0;
    m_E = 0.0;
    m_F = 0.0;
    m_ionDegree = 0.0;
}

void PlasmaPhase::setTemperature(const double temp)
{
    Phase::setTemperature(temp);
    m_kT = Boltzmann * temp / ElectronCharge;
}

void PlasmaPhase::updateElectronEnergyDistribution()
{
    if (m_distributionType == "discretized") {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Invalid for discretized electron energy distribution.");
    } else if (m_distributionType == "isotropic") {
        setIsotropicElectronEnergyDistribution();
    } else if (m_distributionType == "TwoTermApproximation") {
        writelog("call to calculateDistributionFunction()\n");
        auto ierr = ptrEEDFSolver->calculateDistributionFunction();
        if (ierr == 0) {
            auto x = ptrEEDFSolver->getGridEdge();
            auto y = ptrEEDFSolver->getEEDFEdge();
            m_nPoints = x.size();
            m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(x.data(), m_nPoints);
            m_electronEnergyDist = Eigen::Map<const Eigen::ArrayXd>(y.data(), m_nPoints);
        } else {
            throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
                "Call to calculateDistributionFunction failed.");
        }
    }
    electronEnergyDistributionChanged();
    updateElectronTemperatureFromEnergyDist();
    writelog("Done!\n");
}

void PlasmaPhase::normalizeElectronEnergyDistribution() {
    Eigen::ArrayXd eps32 = m_electronEnergyLevels.pow(3./2.);
    double norm = 2./3. * numericalQuadrature(m_quadratureMethod,
                                              m_electronEnergyDist, eps32);
    if (norm < 0.0) {
        throw CanteraError("PlasmaPhase::normalizeElectronEnergyDistribution",
                           "The norm is negative. This might be caused by bad "
                           "electron energy distribution");
    }
    m_electronEnergyDist /= norm;
}

void PlasmaPhase::setElectronEnergyDistributionType(const string& type)
{
    if (type == "discretized" ||
        type == "isotropic" ||
        type == "TwoTermApproximation") {
        m_distributionType = type;
    } else {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistributionType",
            "Unknown type for electron energy distribution.");
    }
}

void PlasmaPhase::setIsotropicElectronEnergyDistribution()
{
    m_electronEnergyDist.resize(m_nPoints);
    double x = m_isotropicShapeFactor;
    double gamma1 = boost::math::tgamma(3.0 / 2.0 * x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 * x);
    double c1 = x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = x * std::pow(gamma2 / gamma1, x);
    m_electronEnergyDist =
        c1 * m_electronEnergyLevels.sqrt() /
        std::pow(meanElectronEnergy(), 1.5) *
        (-c2 * (m_electronEnergyLevels /
        meanElectronEnergy()).pow(x)).exp();
    checkElectronEnergyDistribution();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    m_electronTemp = Te;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setMeanElectronEnergy(double energy) {
    m_electronTemp = 2.0 / 3.0 * energy * ElectronCharge / Boltzmann;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setElectronEnergyLevels(const double* levels, size_t length)
{
    m_nPoints = length;
    m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(levels, length);
    checkElectronEnergyLevels();
    electronEnergyLevelChanged();
    updateElectronEnergyDistribution();
}

void PlasmaPhase::electronEnergyDistributionChanged()
{
    m_distNum++;
}

void PlasmaPhase::electronEnergyLevelChanged()
{
    // The cross sections are interpolated on the energy levels
    if (m_collisions.size() > 0) {
        updateInterpolatedCrossSections();
    }
    m_levelNum++;
}

void PlasmaPhase::checkElectronEnergyLevels() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if (m_electronEnergyLevels[0] < 0.0 || (h <= 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
            "Values of electron energy levels need to be positive and "
            "monotonically increasing.");
    }
}

void PlasmaPhase::checkElectronEnergyDistribution() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if ((m_electronEnergyDist < 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
            "Values of electron energy distribution cannot be negative.");
    }
    if (m_electronEnergyDist[m_nPoints - 1] > 0.01) {
        warn_user("PlasmaPhase::checkElectronEnergyDistribution",
        "The value of the last element of electron energy distribution exceed 0.01. "
        "This indicates that the value of electron energy level is not high enough "
        "to contain the isotropic distribution at mean electron energy of "
        "{} eV", meanElectronEnergy());
    }
}

void PlasmaPhase::setDiscretizedElectronEnergyDist(const double* levels,
                                                const double* dist,
                                                size_t length)
{
    m_distributionType = "discretized";
    m_nPoints = length;
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels, length);
    m_electronEnergyDist =
        Eigen::Map<const Eigen::ArrayXd>(dist, length);
    checkElectronEnergyLevels();
    if (m_do_normalizeElectronEnergyDist) {
        normalizeElectronEnergyDistribution();
    }
    checkElectronEnergyDistribution();
    updateElectronTemperatureFromEnergyDist();
    electronEnergyLevelChanged();
    electronEnergyDistributionChanged();
}

void PlasmaPhase::setDiscretizedElectronEnergyDist(const double* dist,
                                                   size_t length)
{
    m_distributionType = "discretized";
    m_nPoints = length;
    m_electronEnergyDist =
        Eigen::Map<const Eigen::ArrayXd>(dist, length);
    checkElectronEnergyLevels();
    if (m_do_normalizeElectronEnergyDist) {
        normalizeElectronEnergyDistribution();
    }
    checkElectronEnergyDistribution();
    updateElectronTemperatureFromEnergyDist();
    electronEnergyDistributionChanged();
}

void PlasmaPhase::updateElectronTemperatureFromEnergyDist()
{
    // calculate mean electron energy and electron temperature
    Eigen::ArrayXd eps52 = m_electronEnergyLevels.pow(5./2.);
    double epsilon_m = 2.0 / 5.0 * numericalQuadrature(m_quadratureMethod,
                                                       m_electronEnergyDist, eps52);
    m_electronTemp = 2.0 / 3.0 * epsilon_m * ElectronCharge / Boltzmann;
}

void PlasmaPhase::setIsotropicShapeFactor(double x) {
    m_isotropicShapeFactor = x;
    setIsotropicElectronEnergyDistribution();
}

void PlasmaPhase::getParameters(AnyMap& phaseNode) const
{
    IdealGasPhase::getParameters(phaseNode);
    AnyMap eedf;
    eedf["type"] = m_distributionType;
    vector<double> levels(m_nPoints);
    Eigen::Map<Eigen::ArrayXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
    eedf["energy-levels"] = levels;
    if (m_distributionType == "isotropic") {
        eedf["shape-factor"] = m_isotropicShapeFactor;
        eedf["mean-electron-energy"].setQuantity(meanElectronEnergy(), "eV");
    } else if (m_distributionType == "discretized") {
        vector<double> dist(m_nPoints);
        Eigen::Map<Eigen::ArrayXd>(dist.data(), m_nPoints) = m_electronEnergyDist;
        eedf["distribution"] = dist;
        eedf["normalize"] = m_do_normalizeElectronEnergyDist;
    }
    phaseNode["electron-energy-distribution"] = std::move(eedf);
}

void PlasmaPhase::setParameters(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    IdealGasPhase::setParameters(phaseNode, rootNode);
    m_root = rootNode;
    if (phaseNode.hasKey("electron-energy-distribution")) {
        const AnyMap eedf = phaseNode["electron-energy-distribution"].as<AnyMap>();
        m_distributionType = eedf["type"].asString();
        if (m_distributionType == "isotropic") {
            if (eedf.hasKey("shape-factor")) {
                setIsotropicShapeFactor(eedf["shape-factor"].asDouble());
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires shape-factor key.");
            }
            if (eedf.hasKey("mean-electron-energy")) {
                double energy = eedf.convert("mean-electron-energy", "eV");
                setMeanElectronEnergy(energy);
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires electron-temperature key.");
            }
            if (eedf.hasKey("energy-levels")) {
                setElectronEnergyLevels(eedf["energy-levels"].asVector<double>().data(),
                                        eedf["energy-levels"].asVector<double>().size());
            }
            setIsotropicElectronEnergyDistribution();
        } else if (m_distributionType == "discretized") {
            if (!eedf.hasKey("energy-levels")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key energy-levels.");
            }
            if (!eedf.hasKey("distribution")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key distribution.");
            }
            if (eedf.hasKey("normalize")) {
                enableNormalizeElectronEnergyDist(eedf["normalize"].asBool());
            }
            setDiscretizedElectronEnergyDist(eedf["energy-levels"].asVector<double>().data(),
                                             eedf["distribution"].asVector<double>().data(),
                                             eedf["energy-levels"].asVector<double>().size());
        } else if (m_distributionType == "TwoTermApproximation") {
            if (rootNode.hasKey("cross-sections")) {
                // CQM debug
                writelog("I have cross-sections!\n");
                // By default, add all CS from the 'cross-sections' section
                for (const auto& item : rootNode["cross-sections"].asVector<AnyMap>()) {
                    addElectronCrossSection( newElectronCrossSection(item) );
                }
                writelog("m_ncs = {:3d}\n", m_ncs);
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cross section data are required.");
            }
            // CQM use the energy-levels input as initial grid??
            //m_nPoints = eedf["energy-levels"].asVector<double>().size();
            //auto levels = eedf["energy-levels"].asVector<double>().data();
            //m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(levels, m_nPoints);
            ptrEEDFSolver = make_unique<EEDFTwoTermApproximation>(*this);
            // CQM hard coded for now
            // TODO set kTe_max and ncell from user 
            double kTe_max = 60;
            size_t nGridCells = 301;
            // CQM WARNING m_nPoints is nEdges in PlasmaPhase (i.e. nCells+1)
            // It would be nice to change nPoints to nGridEdges
            m_nPoints = nGridCells + 1;
            ptrEEDFSolver->setLinearGrid(kTe_max, nGridCells);
        }
    }
}

bool PlasmaPhase::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    // ecs->validate();
    m_ecss.push_back(ecs);

    m_energyLevels.push_back(ecs->energyLevel);
    m_crossSections.push_back(ecs->crossSection);

    // shift factor
    if (ecs->kind == "ionization") {
        m_shiftFactor.push_back(2);
    } else {
        m_shiftFactor.push_back(1);
    }

    // scattering-in factor
    if (ecs->kind == "ionization") {
        m_inFactor.push_back(2);
    } else if (ecs->kind == "attachment") {
        m_inFactor.push_back(0);
    } else {
        m_inFactor.push_back(1);
    }

    if (ecs->kind == "effective" || ecs->kind == "elastic") {
        for (size_t k = 0; k < m_ncs; k++) {
            if (target(k) == ecs->target)
                if (kind(k) == "elastic" || kind(k) == "effective") {
                    throw CanteraError("PlasmaPhase::addElectronCrossSection"
                                       "Already contains a data of effective/ELASTIC cross section for '{}'.",
                                       ecs->target);
            }
        }
        m_kElastic.push_back(m_ncs);
    } else {
        m_kInelastic.push_back(m_ncs);
    }

    // add one to number of cross sections
    m_ncs++;

    m_f0_ok = false;

    return true;
}

bool PlasmaPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = IdealGasPhase::addSpecies(spec);
    size_t k = m_kk - 1;

    if (spec->composition.find("E") != spec->composition.end() &&
        spec->composition.size() == 1 &&
        spec->composition["E"] == 1) {
        if (m_electronSpeciesIndex == npos) {
            m_electronSpeciesIndex = k;
        } else {
            throw CanteraError("PlasmaPhase::addSpecies",
                               "Cannot add species, {}. "
                               "Only one electron species is allowed.", spec->name);
        }
    }
    return added;
}

void PlasmaPhase::initThermo()
{
    IdealGasPhase::initThermo();
    // check electron species
    if (m_electronSpeciesIndex == npos) {
        throw CanteraError("PlasmaPhase::initThermo",
                           "No electron species found.");
    }

    m_kinetics = newKinetics("bulk");
    m_kinetics->addThermo(shared_from_this());

    vector<shared_ptr<Reaction>> reactions;
    for (AnyMap R : reactionsAnyMapList(*m_kinetics, m_input, m_root)) {
        shared_ptr<Reaction> reaction = newReaction(R, *m_kinetics);

        // Check if the reaction is related to an existing cross-section loaded in the EEDF solver
        if (reaction->type() == "electron-collision-plasma")
        {
            auto rate = std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(reaction->rate());
            for (size_t k = 0; k < m_ncs; k++)
            {
                if (rate->target() == m_ecss[k]->target && rate->product() == m_ecss[k]->product && rate->kind() == m_ecss[k]->kind)
                {
                    rate->set_crossSections(m_ecss[k]->crossSection);
                    rate->set_energyLevels(m_ecss[k]->energyLevel);
                    rate->set_threshold(m_ecss[k]->threshold);
                    rate->set_cs_ok();
                }
            }
            // Check if cross-sections have been defined for the current reaction:
            // Warning: if the cs is not found and does not already exist CRASH
            if (!(rate->get_cs_ok()))
            {   
                writelog("{} -> {}", rate->target(), rate->product());
                throw CanteraError("PlasmaPhase::initThermo",
                                   "Energy levels and cross section are undefined");
            }
        }
        reactions.push_back(reaction);
    }

    // add reactions to kinetics object
    addReactions(*m_kinetics, reactions);

    // init m_collisions and look for elastic collisions
    m_collisions.resize(0);
    size_t i = 0;
    for (shared_ptr<Reaction> reaction : reactions) {
        if (reaction->type() == "electron-collision-plasma") {
            m_collisions.push_back(reaction);

            if (reaction->reactants == reaction->products) {
                // store the indices of elastic collisions
                m_elasticCollisionIndices.push_back(i);
            }
            // only count the electron collision plasma type
            i++;
        }
    }
    updateInterpolatedCrossSections();
}

void PlasmaPhase::updateInterpolatedCrossSections()
{
    for (shared_ptr<Reaction> collision : m_collisions) {
        auto rate = boost::polymorphic_pointer_downcast
            <ElectronCollisionPlasmaRate>(collision->rate());
        vector<double> cs_interp;
        for (double level : m_electronEnergyLevels) {
            cs_interp.push_back(linearInterp(level,
                rate->energyLevels(), rate->crossSections()));
        }
        // Set the interpolated cross section
        rate->setCrossSectionInterpolated(cs_interp);
    }
}

size_t PlasmaPhase::targetSpeciesIndex(shared_ptr<Reaction> R)
{
    if (R->type() != "electron-collision-plasma") {
        throw CanteraError("PlasmaPhase::targetSpeciesIndex",
            "Invalid reaction type. Type electron-collision-plasma is needed.");
    }
    for (const auto& [name, stoich] : R->reactants) {
        if (name != electronSpeciesName()) {
            return speciesIndex(name);
        }
    }
    throw CanteraError("PlasmaPhase::targetSpeciesIndex",
        "No target found. Target cannot be electron.");
}

vector<double> PlasmaPhase::crossSection(shared_ptr<Reaction> reaction)
{
    if (reaction->type() != "electron-collision-plasma") {
        throw CanteraError("PlasmaPhase::crossSection",
            "Invalid reaction type. Type electron-collision-plasma is needed.");
    } else {
        auto rate = boost::polymorphic_pointer_downcast
            <ElectronCollisionPlasmaRate>(reaction->rate());
        std::vector<double> cs_interp;
        for (double level : m_electronEnergyLevels) {
            cs_interp.push_back(linearInterp(level,
                                rate->energyLevels(),
                                rate->crossSections()));
        }
        return cs_interp;
    }
}

double PlasmaPhase::normalizedElasticElectronEnergyLossRate()
{
    double rate = 0.0;
    // calculate dF/dε (forward difference)
    Eigen::ArrayXd dF(nElectronEnergyLevels());
    for (size_t i = 0; i < nElectronEnergyLevels() - 1; i++) {
        dF[i] = (m_electronEnergyDist[i+1] - m_electronEnergyDist[i]) /
                (m_electronEnergyLevels[i+1] - m_electronEnergyLevels[i]);
    }
    dF[nElectronEnergyLevels()-1] = dF[nElectronEnergyLevels()-2];

    for (size_t i : m_elasticCollisionIndices) {
        size_t k = targetSpeciesIndex(m_collisions[i]);
        // get the interpolated cross sections
        auto collision = boost::polymorphic_pointer_downcast
            <ElectronCollisionPlasmaRate>(m_collisions[i]->rate());
        // Map cross sections to Eigen::ArrayXd
        auto cs_array = Eigen::Map<const Eigen::ArrayXd>(
            collision->crossSectionInterpolated().data(),
            collision->crossSectionInterpolated().size()
        );

        double mass_ratio = ElectronMass / molecularWeight(k) * Avogadro;
        rate += mass_ratio * Avogadro * concentration(k) * (
            simpson(1.0 / 3.0 * m_electronEnergyDist.cwiseProduct(
                    cs_array), m_electronEnergyLevels.pow(3.0)) +
            simpson(Boltzmann * temperature() / ElectronCharge *
                    cs_array.cwiseProduct(dF), m_electronEnergyLevels));
    }
    double gamma = sqrt(2 * ElectronCharge / ElectronMass);

    return 2.0 * gamma * rate;
}

void PlasmaPhase::updateThermo() const
{
    IdealGasPhase::updateThermo();
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tempNow = temperature();
    double electronTempNow = electronTemperature();
    size_t k = m_electronSpeciesIndex;
    // If the electron temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tempNow || cached.state2 != electronTempNow) {
        m_spthermo.update_single(k, electronTemperature(),
                &m_cp0_R[k], &m_h0_RT[k], &m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;
    }
    // update the species Gibbs functions
    m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
    // update the nDensity array
    compute_nDensity();
}

void PlasmaPhase::compute_nDensity() const {
    m_nDensity.resize(m_kk);
    const double* Y = massFractions();
    const vector<double>& mw = molecularWeights();
    double dens = density();
    for (size_t j = 0; j < m_kk; j++) {
        m_nDensity[j] = Y[j] / mw[j] * Avogadro * dens;
    }
}

void PlasmaPhase::compute_electronMobility() const {
    if (m_distributionType == "TwoTermApproximation") {
        m_electronMobility = ptrEEDFSolver->getElectronMobility();
    } else {
        throw NotImplementedError("PlasmaPhase::compute_electronMobility");
    }
}

void PlasmaPhase::getVibrationalEnergies(double* const evib) const
{
    copy(m_evib.begin(), m_evib.end(), evib);
}

void PlasmaPhase::setVibrationalEnergies(const double* const evib)
{
    copy(evib, evib + m_nspevib, m_evib.begin());
}

double PlasmaPhase::enthalpy_mole() const {
    double value = IdealGasPhase::enthalpy_mole();
    value += GasConstant * (electronTemperature() - temperature()) *
             moleFraction(m_electronSpeciesIndex) *
             m_h0_RT[m_electronSpeciesIndex];
    return value;
}

void PlasmaPhase::getGibbs_ref(double* g) const
{
    IdealGasPhase::getGibbs_ref(g);
    g[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getStandardVolumes_ref(double* vol) const
{
    IdealGasPhase::getStandardVolumes_ref(vol);
    vol[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEnthalpies(double* hbar) const
{
    IdealGasPhase::getPartialMolarEnthalpies(hbar);
    hbar[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEntropies(double* sbar) const
{
    IdealGasPhase::getPartialMolarEntropies(sbar);
    double logp = log(pressure());
    double logpe = log(electronPressure());
    sbar[m_electronSpeciesIndex] += GasConstant * (logp - logpe);
}

void PlasmaPhase::getPartialMolarIntEnergies(double* ubar) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
    size_t k = m_electronSpeciesIndex;
    ubar[k] = RTe() * (_h[k] - 1.0);
}

void PlasmaPhase::getChemPotentials(double* mu) const
{
    IdealGasPhase::getChemPotentials(mu);
    size_t k = m_electronSpeciesIndex;
    double xx = std::max(SmallNumber, moleFraction(k));
    mu[k] += (RTe() - RT()) * log(xx);
}

void PlasmaPhase::getStandardChemPotentials(double* muStar) const
{
    IdealGasPhase::getStandardChemPotentials(muStar);
    size_t k = m_electronSpeciesIndex;
    muStar[k] -= log(pressure() / refPressure()) * RT();
    muStar[k] += log(electronPressure() / refPressure()) * RTe();
}

void PlasmaPhase::getEntropy_R(double* sr) const
{
    const vector<double>& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            sr[k] -= tmp;
        } else {
            sr[k] -= log(electronPressure() / refPressure());
        }
    }
}

void PlasmaPhase::getGibbs_RT(double* grt) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            grt[k] += tmp;
        } else {
            grt[k] += log(electronPressure() / refPressure());
        }
    }
}

}
