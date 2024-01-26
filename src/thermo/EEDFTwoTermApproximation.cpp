/**
 *  @file EEDFTwoTermApproximation.cpp
 *  EEDF Two-Term approximation solver.  Implementation file for class
 *  EEDFTwoTermApproximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

EEDFTwoTermApproximation::EEDFTwoTermApproximation(PlasmaPhase& s)
{
    initialize(s);
}

void EEDFTwoTermApproximation::initialize(PlasmaPhase& s)
{
    // store a pointer to s.
    m_phase = &s;
}

int EEDFTwoTermApproximation::compute(PlasmaPhase& s, double EN, int loglevel)
{
    // TODO
    // -> call to converge to get f0
    // -> update EEDF dist and grid in the PlasmaPhase object!
    throw CanteraError("EEDFTwoTermApproximation::compute", "To be implemented");
    return 0;
}

void EEDFTwoTermApproximation::converge(Eigen::VectorXd& f0)
{
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
        // writelog("After iteration {:3d}, err = {:.3e} (target: {:.3e}), delta = {:.3e}\n", 
        //                   n + 1, err1, m_rtol, delta);
        // writelog("err1 = {:14.5g} \n",err1);
        if (err1 < options.m_rtol) {
            // writelog("Boltzmann solver convergence after {:d} iterations\n", n);
            break;
        } else if (n == options.m_maxn - 1) {
            throw CanteraError("WeaklyIonizedGas::converge", "Convergence failed");
        }
    }
}

Eigen::VectorXd EEDFTwoTermApproximation::iterate(const Eigen::VectorXd& f0, double delta)
{
    SparseMat_fp PQ(options.m_points, options.m_points);
    vector_fp g = vector_g(f0);
    for (size_t k : m_phase->kInelastic()) {
        PQ += (matrix_Q(g, k) - matrix_P(g, k)) * m_phase->X_targets()[m_phase->klocTargets()[k]];
    }

    SparseMat_fp A = matrix_A(f0);
    SparseMat_fp I(options.m_points, options.m_points);
    for (size_t i = 0; i < options.m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    // Matrix decomposition

    // SimplicialLDLT :
    // Eigen::SimplicialLDLT<SparseMat_fp> solver(A);

    // SparseLU :
    Eigen::SparseLU<SparseMat_fp> solver(A);
    
    // SparseQR :
    // Eigen::SparseQR<SparseMat_fp, Eigen::COLAMDOrdering<int>> solver;
    // Eigen::SparseQR<SparseMat_fp, Eigen::AMDOrdering<int>> solver;
    // solver.compute(A);

    if (solver.info() != Eigen::Success) {
        //writelog("Decomposition failed\n");
        throw CanteraError("EEDFTwoTermApproximation::iterate", "Decomposition failed");
        return f0;
    }

    // solve f0
    Eigen::VectorXd f1 = solver.solve(f0);
    if(solver.info() != Eigen::Success) {
        //writelog("Solving failed \n");
        throw CanteraError("EEDFTwoTermApproximation::iterate", "Solving failed");
        return f0;
    }

    f1 /= norm(f1, m_gridCenter);
    return f1;
}

double EEDFTwoTermApproximation::integralPQ(double a, double b, double u0, double u1,
                                            double g, double x0)
{
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
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double j = m_j[k][n];
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
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double i = m_i[k][n];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double q = m_inFactor[k] * m_gamma * r;
        tripletList.push_back(Triplet_fp(i, j, q));
    }
    SparseMat_fp Q(options.m_points, options.m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_A(const Eigen::VectorXd& f0)
{
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
    }

    double alpha;
    if (options.m_growth == "spatial") {
        double mu = electronMobility(m_f0);
        double D = electronDiffusivity(m_f0);
        alpha = (mu * m_phase->E() - sqrt(pow(mu * m_phase->E(), 2) - 4 * D * nu * m_phase->N())) / 2.0 / D / m_phase->N();
    }


    for (size_t j = 1; j < options.m_points; j++) {
        double sigma_tilde;
        if (options.m_growth == "temporal") {
            sigma_tilde = m_totalCrossSectionEdge[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
        }
        else {
            sigma_tilde = m_totalCrossSectionEdge[j];
        }
        double omega = 2 * Pi * m_phase->F();
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

    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        if (m_phase->kind(k) == "ionization" ||
            m_phase->kind(k) == "attachment") {
            SparseMat_fp PQ = (matrix_Q(g, k) - matrix_P(g, k)) *
                              m_phase->X_targets()[m_phase->klocTargets()[k]];
            Eigen::VectorXd s = PQ * f0;
            for (size_t i = 0; i < options.m_points; i++) {
                nu += s[i];
            }
        }
    }
    return nu;
}

double EEDFTwoTermApproximation::electronDiffusivity(const Eigen::VectorXd m_f0)
{
    vector_fp y(options.m_points, 0.0);
    double nu = netProductionFreq(m_f0);
    for (size_t i = 0; i < options.m_points; i++) {
        if (m_gridCenter[i] != 0.0) {
            y[i] = m_gridCenter[i] * m_f0(i) /
                   (m_totalCrossSectionCenter[i] + nu / m_gamma / pow(m_gridCenter[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridCenter.data(), m_gridCenter.size());
    return 1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

double EEDFTwoTermApproximation::electronMobility(const Eigen::VectorXd m_f0)
{
    double nu = netProductionFreq(m_f0);
    vector_fp y(options.m_points + 1, 0.0);
    for (size_t i = 1; i < options.m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            y[i] = m_gridEdge[i] * df0 /
                   (m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridEdge.data(), m_gridEdge.size());
    return -1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

double EEDFTwoTermApproximation::norm(const Eigen::VectorXd& f, const vector_fp& grid)
{
    // TODO
    // CQM /!\ already a norm function in PlasmaPhase
    throw CanteraError("EEDFTwoTermApproximation::norm", "To be implemented");
    return 0.0;
}

} // end of namespace Cantera