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
            delta *= std::log(options.m_factorM) / (std::log(err0) - std::log(err1));
        }
        Eigen::VectorXd f0_old = f0;
        f0 = iterate(f0_old, delta);
        err0 = err1;
        Eigen::VectorXd Df0(options.m_points);
        for (size_t i = 0; i < options.m_points; i++) {
            Df0(i) = std::abs(f0_old(i) - f0(i));
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

vector_fp EEDFTwoTermApproximation::vector_g(const Eigen::VectorXd& f0)
{
    // TODO
    vector_fp g(options.m_points, 0.0);
    return g;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_P(const vector_fp& g, size_t k)
{
    // TODO
    SparseMat_fp P(options.m_points, options.m_points);
    return P;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_Q(const vector_fp& g, size_t k)
{
    // TODO
    SparseMat_fp Q(options.m_points, options.m_points);
    return Q;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_A(const Eigen::VectorXd& f0)
{
    // TODO
    SparseMat_fp A(options.m_points, options.m_points);
    SparseMat_fp G(options.m_points, options.m_points);
    return A + G;
}

double EEDFTwoTermApproximation::norm(const Eigen::VectorXd& f, const vector_fp& grid)
{
    // TODO
    // CQM /!\ already a norn function in PlasmaPhase
    throw CanteraError("EEDFTwoTermApproximation::norm", "To be implemented");
    return 0.0;
}

} // end of namespace Cantera