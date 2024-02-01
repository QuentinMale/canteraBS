/**
 *  @file EEDFTwoTermApproximation.cpp
 *  EEDF Two-Term approximation solver.  Implementation file for class
 *  EEDFTwoTermApproximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/PlasmaPhase.h"
#include <iostream>

namespace Cantera
{

EEDFTwoTermApproximation::EEDFTwoTermApproximation(PlasmaPhase& s)
{
    writelog("EEDFTwoTermApproximation\n");
    initialize(s);
}

void EEDFTwoTermApproximation::initialize(PlasmaPhase& s)
{
    writelog("initialize EEDFTwoTermApproximation\n");
    // store a pointer to s.
    m_phase = &s;
    m_first_call = true;
    m_has_EEDF = false;
    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

void EEDFTwoTermApproximation::setLinearGrid(double& kTe_max, size_t& ncell)
{
    writelog("Grid info : Linear grid is used \n");
    writelog("Grid info : Maximum energy of the grid is {:15.3g} [eV]\n", kTe_max);
    options.m_points = ncell;
    m_gridCenter.resize(options.m_points);
    m_gridEdge.resize(options.m_points + 1);
    m_f0.resize(options.m_points);
    m_f0_edge.resize(options.m_points + 1);
    for (size_t j = 0; j < options.m_points; j++) {
        m_gridCenter[j] = kTe_max * ( j + 0.5 ) / options.m_points;
        m_gridEdge[j] = kTe_max * j / options.m_points;
    }
    m_gridEdge[options.m_points] = kTe_max;
    setGridCache();
}

int EEDFTwoTermApproximation::calculateDistributionFunction()
{
    // TODO
    // -> call to converge to get f0
    // -> update EEDF dist and grid in the PlasmaPhase object!

    // During the first call to this function the indices of target species need to be defined
    if (m_first_call)
    {
        writelog("First call to calculateDistributionFunction\n");
        initSpeciesIndexCS();
        m_first_call = false;
    } else {
        writelog("pass init\n");
    }

    update_mole_fractions();
    checkSpeciesNoCrossSection();
    updateCS();

    if (!m_has_EEDF) {
        if (options.m_firstguess == "maxwell") {
            writelog("First guess EEDF maxwell\n");
            //auto kTe_max = 10.0 * options.m_init_kTe;
            for (size_t j = 0; j < options.m_points; j++) {
                m_f0(j) = 2.0 * pow(1.0 / Pi, 0.5) * pow(options.m_init_kTe, -3. / 2.) *
                          exp(-m_gridCenter[j] / options.m_init_kTe);
            }
        } else {
            throw CanteraError("EEDFTwoTermApproximation::calculateDistributionFunction",
                               " unknown EEDF first guess");
        }
    }

    // Start of monitoring
    m_timer_eedf->start();

    // Computation of the EEDF
    converge(m_f0);

    // End of monitoring
    m_timer_eedf->stop();

    // write the EEDF at grid edges
    vector<double> f(m_f0.data(), m_f0.data() + m_f0.rows() * m_f0.cols());
    vector<double> x(m_gridCenter.data(), m_gridCenter.data() + m_gridCenter.rows() * m_gridCenter.cols());
    for (size_t i = 0; i < options.m_points + 1; i++) {
        m_f0_edge[i] = linearInterp(m_gridEdge[i], x, f);
    }

    m_has_EEDF = true;

    return 0;

}

void EEDFTwoTermApproximation::converge(Eigen::VectorXd& f0)
{
    writelog("EEDFTwoTermApproximation::converge\n");
    double err0 = 0.0;
    double err1 = 0.0;
    double delta = options.m_delta0;
    for (size_t n = 0; n < options.m_maxn; n++) {
        if (0.0 < err1 && err1 < err0) {
            // log extrapolation attempting to reduce the error a factor m
            delta *= log(options.m_factorM) / (log(err0) - log(err1));
        }
        Eigen::VectorXd f0_old = f0;
        f0 = iterate(f0_old, delta);
        err0 = err1;
        Eigen::VectorXd Df0(options.m_points);
        for (size_t i = 0; i < options.m_points; i++) {
            Df0(i) = abs(f0_old(i) - f0(i));
        }
        err1 = norm(Df0, m_gridCenter);
        writelog("After iteration {:3d}, err = {:.3e} (target: {:.3e}), delta = {:.3e}\n", 
                  n + 1, err1, options.m_rtol, delta);
        writelog("err1 = {:14.5g} \n",err1);
        if (err1 < options.m_rtol) {
            writelog("Boltzmann solver convergence after {:d} iterations\n", n);
            break;
        } else if (n == options.m_maxn - 1) {
            throw CanteraError("WeaklyIonizedGas::converge", "Convergence failed");
        }
    }
}

Eigen::VectorXd EEDFTwoTermApproximation::iterate(const Eigen::VectorXd& f0, double delta)
{
    // CQM multiple call to vector_* and matrix_*
    // probably extremely ineficient
    // must be refactored!!

    writelog("EEDFTwoTermApproximation::iterate\n");
    SparseMat_fp PQ(options.m_points, options.m_points);
    vector_fp g = vector_g(f0);
    for (size_t k : m_phase->kInelastic()) {
        PQ += (matrix_Q(g, k) - matrix_P(g, k)) * m_X_targets[m_klocTargets[k]];
    }

    SparseMat_fp A = matrix_A(f0);
    SparseMat_fp I(options.m_points, options.m_points);
    for (size_t i = 0; i < options.m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    // Check matrix validity
    //writelog("{:d}rows {:d}cols\n", A.rows(), A.cols());
    writelog("Number of non zero values: {:d}\n", A.nonZeros());
    if (!A.isVector()) {
        writelog("The matrix A is not sparse!\n");
    }
    if (!A.isCompressed()) {
        writelog("The matrix A is not in compressed form!\n");
    }

    // Matrix decomposition

    // SparseLU :
    Eigen::SparseLU<SparseMat_fp> solver(A);
    if (solver.info() == Eigen::NumericalIssue) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver: NumericalIssue");
    } else if (solver.info() == Eigen::InvalidInput) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver: InvalidInput");
    }
    if (solver.info() != Eigen::Success) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver", "Decomposition failed");
        return f0;
    }

    // solve f0
    Eigen::VectorXd f1 = solver.solve(f0);
    if(solver.info() != Eigen::Success) {
        throw CanteraError("EEDFTwoTermApproximation::iterate", "Solving failed");
        return f0;
    }

    f1 /= norm(f1, m_gridCenter);
    return f1;
}

double EEDFTwoTermApproximation::integralPQ(double a, double b, double u0, double u1,
                                            double g, double x0)
{
    //writelog("EEDFTwoTermApproximation::integralPQ\n");
    double A1;
    double A2;
    if (g != 0.0) {
        double expm1a = expm1(g * (-a + x0));
        double expm1b = expm1(g * (-b + x0));
        double ag = a * g;
        double ag1 = ag + 1;
        double bg = b * g;
        double bg1 = bg + 1;
        A1 = (expm1a * ag1 + ag - expm1b * bg1 - bg) / (g*g);
        A2 = (expm1a * (2 * ag1 + ag * ag) + ag * (ag + 2) -
              expm1b * (2 * bg1 + bg * bg) - bg * (bg + 2)) / (g*g*g);
    } else {
        A1 = 0.5 * (b*b - a*a);
        A2 = 1.0 / 3.0 * (b*b*b - a*a*a);
    }

    // The interpolation formula of u(x) = c0 + c1 * x
    double c0 = (a * u1 - b * u0) / (a - b);
    double c1 = (u0 - u1) / (a - b);

    return c0 * A1 + c1 * A2;
}

vector_fp EEDFTwoTermApproximation::vector_g(const Eigen::VectorXd& f0)
{
    //writelog("EEDFTwoTermApproximation::vector_g\n");
    vector_fp g(options.m_points, 0.0);
    g[0] = log(f0(1)/f0(0)) / (m_gridCenter[1] - m_gridCenter[0]);
    size_t N = options.m_points - 1;
    g[N] = log(f0(N)/f0(N-1)) / (m_gridCenter[N] - m_gridCenter[N-1]);
    for (size_t i = 1; i < options.m_points - 1; i++) {
        g[i] = log(f0(i+1)/f0(i-1)) / (m_gridCenter[i+1] - m_gridCenter[i-1]);
    }
    return g;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_P(const vector_fp& g, size_t k)
{
    //writelog("EEDFTwoTermApproximation::matrix_P\n");
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        size_t j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double p = m_gamma * r;
        tripletList.push_back(Triplet_fp(j, j, p));
    }
    SparseMat_fp P(options.m_points, options.m_points);
    P.setFromTriplets(tripletList.begin(), tripletList.end());
    return P;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_Q(const vector_fp& g, size_t k)
{
    //writelog("EEDFTwoTermApproximation::matrix_Q\n");
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        size_t i = m_i[k][n];
        size_t j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double q = m_phase->inFactor()[k] * m_gamma * r;
        tripletList.push_back(Triplet_fp(i, j, q));
    }
    SparseMat_fp Q(options.m_points, options.m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_A(const Eigen::VectorXd& f0)
{
    //writelog("EEDFTwoTermApproximation::matrix_A\n");
    vector_fp a0(options.m_points + 1);
    vector_fp a1(options.m_points + 1);
    size_t N = options.m_points - 1;
    // Scharfetter-Gummel scheme
    double nu = netProductionFreq(f0);
    a0[0] = NAN;
    a1[0] = NAN;
    a0[N+1] = NAN;
    a1[N+1] = NAN;

    // Electron-electron collisions declarations
    double a;
    vector_fp A1, A2, A3;
    if (m_eeCol) {
        // TODO
        //eeColIntegrals(A1, A2, A3, a, options.m_points);
        throw CanteraError("EEDFTwoTermApproximation::matrix_A",
            "eeCol to be implemented");
    }

    double alpha;
    if (options.m_growth == "spatial") {
        double mu = electronMobility(f0);
        double D = electronDiffusivity(f0);
        alpha = (mu * m_phase->E() - sqrt(pow(mu * m_phase->E(), 2) - 4 * D * nu * m_phase->N())) / 2.0 / D / m_phase->N();
    } else {
        alpha = 0.0;
    }

    double sigma_tilde;
    double omega = 2 * Pi * m_phase->F();
    for (size_t j = 1; j < options.m_points; j++) {
        if (options.m_growth == "temporal") {
            sigma_tilde = m_totalCrossSectionEdge[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
        }
        else {
            sigma_tilde = m_totalCrossSectionEdge[j];
        }
        double q = omega / (m_phase->N() * m_gamma * pow(m_gridEdge[j], 0.5));
        double W = -m_gamma * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double F = sigma_tilde * sigma_tilde / (sigma_tilde * sigma_tilde + q * q);
        double DA = m_gamma / 3.0 * pow(m_phase->E() / m_phase->N(), 2.0) * m_gridEdge[j];
        double DB = m_gamma * m_phase->kT() * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double D = DA / sigma_tilde * F + DB;
        if (m_eeCol) {
            W -= 3 * a * m_phase->ionDegree() * A1[j];
            D += 2 * a * m_phase->ionDegree() * (A2[j] + pow(m_gridEdge[j], 1.5) * A3[j]);
        }
        if (options.m_growth == "spatial") {
            W -= m_gamma / 3.0 * 2 * alpha * m_phase->E() / m_phase->N() * m_gridEdge[j] / sigma_tilde;
        }
        double z = W * (m_gridCenter[j] - m_gridCenter[j-1]) / D;
        a0[j] = W / (1 - std::exp(-z));
        a1[j] = W / (1 - std::exp(z));
    }

    std::vector<Triplet_fp> tripletList;
    // center diagonal
    // zero flux b.c. at energy = 0
    tripletList.push_back(Triplet_fp(0, 0, a0[1]));

    for (size_t j = 1; j < options.m_points - 1; j++) {
        tripletList.push_back(Triplet_fp(j, j, a0[j+1] - a1[j]));
    }

    // upper diagonal
    for (size_t j = 0; j < options.m_points - 1; j++) {
        tripletList.push_back(Triplet_fp(j, j+1, a1[j+1]));
    }

    // lower diagonal
    for (size_t j = 1; j < options.m_points; j++) {
        tripletList.push_back(Triplet_fp(j, j-1, -a0[j]));
    }

    // zero flux b.c.
    tripletList.push_back(Triplet_fp(N, N, -a1[N]));

    SparseMat_fp A(options.m_points, options.m_points);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    //plus G
    SparseMat_fp G(options.m_points, options.m_points);
    if (options.m_growth == "temporal") {
        for (size_t i = 0; i < options.m_points; i++) {
            G.insert(i, i) = 2.0 / 3.0 * (pow(m_gridEdge[i+1], 1.5) - pow(m_gridEdge[i], 1.5)) * nu;
        }
    }
    else if (options.m_growth == "spatial") {
        for (size_t i = 0; i < options.m_points; i++) {
            double sigma_c = 0.5 * (m_totalCrossSectionEdge[i] + m_totalCrossSectionEdge[i + 1]);
            G.insert(i, i) = - alpha * m_gamma / 3 * (alpha * (pow(m_gridEdge[i + 1], 2) - pow(m_gridEdge[i], 2)) / sigma_c / 2
                 - m_phase->E() / m_phase->N() * (m_gridEdge[i + 1] / m_totalCrossSectionEdge[i + 1] - m_gridEdge[i] / m_totalCrossSectionEdge[i]));
        }
    }

    return A + G;
}

double EEDFTwoTermApproximation::netProductionFreq(const Eigen::VectorXd& f0)
{
    double nu = 0.0;
    vector_fp g = vector_g(f0);

    //writelog("EEDFTwoTermApproximation::netProductionFreq\n");

    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        if (m_phase->kind(k) == "ionization" ||
            m_phase->kind(k) == "attachment") {
            SparseMat_fp PQ = (matrix_Q(g, k) - matrix_P(g, k)) *
                              m_X_targets[m_klocTargets[k]];
            Eigen::VectorXd s = PQ * f0;
            for (size_t i = 0; i < options.m_points; i++) {
                nu += s[i];
            }
        }
    }
    return nu;
}

double EEDFTwoTermApproximation::electronDiffusivity(const Eigen::VectorXd& f0)
{
    vector_fp y(options.m_points, 0.0);
    double nu = netProductionFreq(f0);
    for (size_t i = 0; i < options.m_points; i++) {
        if (m_gridCenter[i] != 0.0) {
            y[i] = m_gridCenter[i] * f0(i) /
                   (m_totalCrossSectionCenter[i] + nu / m_gamma / pow(m_gridCenter[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridCenter.data(), m_gridCenter.size());
    return 1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

double EEDFTwoTermApproximation::electronMobility(const Eigen::VectorXd& f0)
{
    double nu = netProductionFreq(f0);
    vector_fp y(options.m_points + 1, 0.0);
    for (size_t i = 1; i < options.m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (f0(i) - f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            y[i] = m_gridEdge[i] * df0 /
                   (m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridEdge.data(), m_gridEdge.size());
    return -1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

void EEDFTwoTermApproximation::initSpeciesIndexCS()
{
    //writelog("initSpeciesIndexCS 1\n");
    // set up target index
    m_kTargets.resize(m_phase->nElectronCrossSections());
    m_klocTargets.resize(m_phase->nElectronCrossSections());
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++)
    {
        writelog("{:d} {:s}\n", k, m_phase->target(k));
        m_kTargets[k] = m_phase->speciesIndex(m_phase->target(k));
        if (m_kTargets[k] == string::npos) {
            throw CanteraError("EEDFTwoTermApproximation::initSpeciesIndexCS"
                               " species not found!",
                               m_phase->target(k));
        }
        // Check if it is a new target or not :
        auto it = find(m_k_lg_Targets.begin(), m_k_lg_Targets.end(), m_kTargets[k]);
        writelog("Indice = {:d}\n", distance(m_k_lg_Targets.begin(), it));
        if (it == m_k_lg_Targets.end()){
            writelog("New target found: {:s} with index {:d}\n", m_phase->target(k), distance(m_k_lg_Targets.begin(), it));
            m_k_lg_Targets.push_back(m_kTargets[k]);
            m_klocTargets[k] = m_k_lg_Targets.size() - 1;
        } else {
            m_klocTargets[k] = distance(m_k_lg_Targets.begin(), it);
        }
    }

    //writelog("initSpeciesIndexCS 2\n");

    writelog("Number of target species found: {:d}\n", m_k_lg_Targets.size());
    m_X_targets.resize(m_k_lg_Targets.size());
    m_X_targets_prev.resize(m_k_lg_Targets.size());
    writelog("Number of target species found: {:d}\n", m_X_targets.size());
    //writelog("initSpeciesIndexCS 21\n");
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        //writelog("initSpeciesIndexCS 22 {:d}\n", k);
        writelog("m_k_lg_Targets[{:d}] = {:d}\n", k, m_k_lg_Targets[k]);
        writelog("moleFraction[{:d}] = {:g}\n", m_k_lg_Targets[k], m_phase->moleFraction(m_k_lg_Targets[k]));
        size_t k_glob = m_k_lg_Targets[k];
        m_X_targets[k] = m_phase->moleFraction(k_glob);
        m_X_targets_prev[k] = m_phase->moleFraction(k_glob);
        writelog("The target number {:d} has X = {:.3g}\n", k, m_X_targets[k]);
    }

    //writelog("initSpeciesIndexCS 3\n");

    // set up indices of species which has no cross-section data
    for (size_t k = 0; k < m_phase->nSpecies(); k++)
    {
        auto it = std::find(m_kTargets.begin(), m_kTargets.end(), k);
        if (it == m_kTargets.end()) {
            m_kOthers.push_back(k);
        }
    }
}

void EEDFTwoTermApproximation::checkSpeciesNoCrossSection()
{
    // warn that a specific species needs cross-section data.
    for (size_t k : m_kOthers) {
        if (m_phase->moleFraction(k) > options.m_moleFractionThreshold) {
            writelog("EEDFTwoTermApproximation:checkSpeciesNoCrossSection\n");
            writelog("Warning:The mole fraction of species {} is more than 0.01 (X = {:.3g}) but it has no cross-section data\n", m_phase->speciesName(k), m_phase->moleFraction(k));
        }
    }
}

void EEDFTwoTermApproximation::updateCS()
{
    writelog("updateCS\n");
    // Compute sigma_m and sigma_\epsilon
    calculateTotalCrossSection();
    calculateTotalElasticCrossSection();
}

// Update the species mole fractions used for EEDF computation
void EEDFTwoTermApproximation::update_mole_fractions()
{
    writelog("Update mole fractions in EEDFTwoTermApproximation\n");
    double tmp_sum = 0.0;
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        writelog("The target number {:d} has X = {:.3g}\n", k, m_phase->moleFraction(m_k_lg_Targets[k]));
        writelog("update X {:d}\n", k);
        m_X_targets[k] = m_phase->moleFraction(m_k_lg_Targets[k]);
        tmp_sum = tmp_sum + m_phase->moleFraction(m_k_lg_Targets[k]);
    }
    //writelog("Update mole fractions in EEDFTwoTermApproximation 2\n");
    writelog("Sum of mole fraction is equal to {:.2g}\n", tmp_sum);

    // Normalize the mole fractions to unity:
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        m_X_targets[k] = m_X_targets[k] / tmp_sum;
        writelog("The target number {:d} has X = {:.3g}\n", k, m_X_targets[k]);
        // if (fabs(m_X_targets[k] - m_X_targets_prev[k]) >= m_X_atol)
        // {
        //     writelog("Mole fractions change a lot, m_f0_ok is set to false\n");
        //     m_f0_ok = false;
        // }
    }
    //writelog("Update mole fractions in EEDFTwoTermApproximation 3\n");
}

void EEDFTwoTermApproximation::calculateTotalCrossSection()
{
    writelog("calculateTotalCrossSection\n");
    m_totalCrossSectionCenter.assign(options.m_points, 0.0);
    m_totalCrossSectionEdge.assign(options.m_points + 1, 0.0);
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        vector_fp x = m_phase->energyLevels()[k];
        vector_fp y = m_phase->crossSections()[k];
        writelog("Kind :    {:s}\n", m_phase->kind(k));
        writelog("Target :    {:s}\n", m_phase->target(k));
        writelog("Product :    {:s}\n", m_phase->product(k));
        writelog("Check X: {:g} =? {:g}\n", m_phase->moleFraction(m_kTargets[k]), m_X_targets[m_klocTargets[k]]);
        writelog("\n");
        for (size_t i = 0; i < options.m_points; i++) {
            m_totalCrossSectionCenter[i] += m_X_targets[m_klocTargets[k]] *
                                            linearInterp(m_gridCenter[i], x, y);
        }
        for (size_t i = 0; i < options.m_points + 1; i++) {
            m_totalCrossSectionEdge[i] += m_X_targets[m_klocTargets[k]] *
                                          linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::calculateTotalElasticCrossSection()
{
    writelog("calculateTotalElasticCrossSection\n");
    m_sigmaElastic.clear();
    m_sigmaElastic.resize(options.m_points, 0.0);
    for (size_t k : m_phase->kElastic()) {
        vector_fp x = m_phase->energyLevels()[k];
        vector_fp y = m_phase->crossSections()[k];
        // Note:
        // moleFraction(m_kTargets[k]) <=> m_X_targets[m_klocTargets[k]]
        double mass_ratio = ElectronMass / (m_phase->molecularWeight(m_kTargets[k]) / Avogadro);
        for (size_t i = 0; i < options.m_points; i++) {
            m_sigmaElastic[i] += 2.0 * mass_ratio * m_X_targets[m_klocTargets[k]] *
                                 linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::setGridCache()
{
    writelog("EEDFTwoTermApproximation::setGridCache\n");
    m_sigma.clear();
    m_sigma.resize(m_phase->nElectronCrossSections());
    m_sigma_offset.clear();
    m_sigma_offset.resize(m_phase->nElectronCrossSections());
    m_eps.clear();
    m_eps.resize(m_phase->nElectronCrossSections());
    m_j.clear();
    m_j.resize(m_phase->nElectronCrossSections());
    m_i.clear();
    m_i.resize(m_phase->nElectronCrossSections());
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        //vector_fp& x = m_phase->energyLevels()[k];
        auto x = m_phase->energyLevels()[k];
        //vector_fp& y = m_phase->crossSections()[k];
        auto y = m_phase->crossSections()[k];
        vector_fp eps1(options.m_points + 1);
        for (size_t i = 0; i < options.m_points + 1; i++) {
            eps1[i] = clip(m_phase->shiftFactor()[k] * m_gridEdge[i] + m_phase->threshold(k),
                           m_gridEdge[0] + 1e-9, m_gridEdge[options.m_points] - 1e-9);
        }
        vector_fp nodes = eps1;
        for (size_t i = 0; i < options.m_points + 1; i++) {
            if (m_gridEdge[i] >= eps1[0] && m_gridEdge[i] <= eps1[options.m_points]) {
                nodes.push_back(m_gridEdge[i]);
            }
        }
        for (size_t i = 0; i < x.size(); i++) {
            if (x[i] >= eps1[0] && x[i] <= eps1[options.m_points]) {
                nodes.push_back(x[i]);
            }
        }

        std::sort(nodes.begin(), nodes.end());
        auto last = std::unique(nodes.begin(), nodes.end());
        nodes.resize(std::distance(nodes.begin(), last));
        vector_fp sigma0(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            sigma0[i] = linearInterp(nodes[i], x, y);
        }

        // search position of cell j
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(m_gridEdge.begin(), m_gridEdge.end(), nodes[i]);
            m_j[k].push_back(low - m_gridEdge.begin() - 1);
        }

        // search position of cell i
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(eps1.begin(), eps1.end(), nodes[i]);
            m_i[k].push_back(low - eps1.begin() - 1);
        }

        // construct sigma
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp sigma{sigma0[i], sigma0[i+1]};
            m_sigma[k].push_back(sigma);
        }

        // construct eps
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp eps{nodes[i], nodes[i+1]};
            m_eps[k].push_back(eps);
        }

        // construct sigma_offset
        auto x_offset = m_phase->energyLevels()[k];
        for (auto& element : x_offset) {
            element -= m_phase->threshold(k);
        }
        for (size_t i = 0; i < options.m_points; i++) {
            m_sigma_offset[k].push_back(linearInterp(m_gridCenter[i], x_offset, y));
        }
    }
}

double EEDFTwoTermApproximation::norm(const Eigen::VectorXd& f, const Eigen::VectorXd& grid)
{
    string m_quadratureMethod = "simpson";
    Eigen::VectorXd p(f.size());
    for (int i = 0; i < f.size(); i++) {
        p[i] = f(i) * pow(grid[i], 0.5);
    }
    return numericalQuadrature(m_quadratureMethod, p, grid);
}

} // end of namespace Cantera
