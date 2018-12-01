//! @file IonGasTransport.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/IonGasTransport.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/stringUtils.h"
#include "MMCollisionInt.h"

namespace Cantera
{
IonGasTransport::IonGasTransport() :
    m_kElectron(npos)
{
    // effective electron collisional cross section to major species
    // The order is N2 O2 CO2 CO H2O H2 CH4
    const std::vector<vector_fp> electron_energy = {
        // N2
        {0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.0085, 0.01, 0.015, 0.02,
        0.03, 0.04, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35,
        0.4, 0.5, 0.7, 1.0, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.5, 2.8, 3.0,
        3.3, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0,
        25.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0,
        1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 10000.0},
        // O2
        {0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.0085, 0.01, 0.015, 0.02,
        0.03, 0.04, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35,
        0.4, 0.5, 0.7, 1.0, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.5, 2.8, 3.0,
        3.3, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0,
        25.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0,
        1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 10000.0},
        // CO2
        {0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.0085, 0.01, 0.015, 0.02,
        0.03, 0.04, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35,
        0.4, 0.5, 0.7, 1.0, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.5, 2.8, 3.0,
        3.3, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0,
        25.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0,
        1000.0},
        // CO
        {0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.0085, 0.01, 0.015, 0.02,
        0.03, 0.04, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35,
        0.4, 0.5, 0.7, 1.0, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.5, 2.8, 3.0,
        3.3, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0,
        25.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0,
        1000.0},
        // H2O
        {0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.14, 0.2, 0.3, 0.5, 0.6, 0.7,
        0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 13.0,
        17.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0, 500.0},
        // H2
        {0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.0085, 0.01, 0.015, 0.02,
        0.03, 0.04, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.35,
        0.4, 0.5, 0.7, 1.0, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.5, 2.8, 3.0,
        3.3, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0,
        25.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0, 700.0,
        1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 10000.0},
        // CH4
        {0.0, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.32, 0.43, 0.55, 0.9, 1.3,
        2.0, 3.0, 6.0, 8.0, 10.0, 20.0, 40.0, 70.0, 100.0, 300.0, 1000.0}
    };
    const std::vector<vector_fp> cross_section = {
        // N2
        {1.1e-20, 1.36e-20, 1.49e-20, 1.62e-20, 1.81e-20, 2e-20, 2.1e-20,
        2.19e-20, 2.55e-20, 2.85e-20, 3.4e-20, 3.85e-20, 4.33e-20, 5.1e-20,
        5.95e-20, 6.45e-20, 7.1e-20, 7.4e-20, 7.9e-20, 8.5e-20, 9e-20, 9.4e-20,
        9.7e-20, 9.9e-20, 1e-19, 1e-19, 1.04e-19, 1.1e-19, 1.2e-19, 1.38e-19,
        1.96e-19, 2.7e-19, 2.85e-19, 3e-19, 2.8e-19, 2.17e-19, 1.72e-19,
        1.47e-19, 1.26e-19, 1.13e-19, 1.09e-19, 1.04e-19, 1.01e-19, 1e-19,
        1.04e-19, 1.09e-19, 1.1e-19, 1.07e-19, 1.02e-19, 9.5e-20, 9e-20,
        8.6e-20, 6.6e-20, 5.8e-20, 4.9e-20, 4.2e-20, 3.3e-20, 2.44e-20,
        1.96e-20, 1.55e-20, 1.12e-20, 8.1e-21, 6.3e-21, 4e-21, 2.9e-21,
        2.1e-21},
        // O2
        {3.5e-21, 3.5e-21, 3.6e-21, 4e-21, 5e-21, 5.8e-21, 6.4e-21, 7e-21,
        8.7e-21, 9.9e-21, 1.24e-20, 1.44e-20, 1.6e-20, 2.1e-20, 2.5e-20,
        2.8e-20, 3.1e-20, 3.3e-20, 3.6e-20, 4.1e-20, 4.5e-20, 4.7e-20, 5.2e-20,
        5.7e-20, 6.1e-20, 7.2e-20, 7.9e-20, 7.9e-20, 7.6e-20, 7.3e-20, 6.9e-20,
        6.6e-20, 6.5e-20, 6.1e-20, 5.8e-20, 5.7e-20, 5.5e-20, 5.45e-20,
        5.5e-20, 5.55e-20, 5.6e-20, 6e-20, 6.6e-20, 7.1e-20, 8e-20, 8.5e-20,
        8.8e-20, 8.7e-20, 8.6e-20, 8.2e-20, 8e-20, 7.7e-20, 6.8e-20, 6.5e-20,
        6.7e-20, 6e-20, 4.9e-20, 3.6e-20, 2.9e-20, 2.12e-20, 1.48e-20,
        1.14e-20, 7.9e-21, 5.1e-21, 3.8e-21, 2.8e-21},
        // CO2
        {6e-18, 5.4e-18, 3.8e-18, 3.07e-18, 2.37e-18, 2e-18, 1.82e-18, 1.7e-18,
        1.38e-18, 1.2e-18, 9.7e-19, 8.5e-19, 7.6e-19, 6.3e-19, 5e-19, 4.4e-19,
        3.9e-19, 3.4e-19, 2.9e-19, 2.4e-19, 1.8e-19, 1.5e-19, 1.3e-19, 1e-19,
        7.1e-20, 5.2e-20, 4.8e-20, 4.7e-20, 4.65e-20, 4.65e-20, 4.85e-20,
        5.05e-20, 5.2e-20, 6.4e-20, 7.6e-20, 9e-20, 1.15e-19, 1.4e-19,
        1.52e-19, 1.48e-19, 1.27e-19, 1e-19, 1e-19, 1.08e-19, 1.21e-19,
        1.31e-19, 1.44e-19, 1.5e-19, 1.58e-19, 1.6e-19, 1.58e-19, 1.26e-19,
        9.5e-20, 8e-20, 6e-20, 4e-20, 3.7e-20, 2.5e-20, 2e-20, 1.6e-20},
        // CO
        {6e-19, 4e-19, 2.5e-19, 1.77e-19, 1.23e-19, 9.8e-20, 8.6e-20, 7.8e-20,
        6.5e-20, 5.9e-20, 5.4e-20, 5.2e-20, 5.4e-20, 6.1e-20, 7.3e-20, 7.7e-20,
        8.8e-20, 9.3e-20, 1e-19, 1.12e-19, 1.21e-19, 1.3e-19, 1.385e-19,
        1.54e-19, 1.65e-19, 1.85e-19, 2.8e-19, 3.7e-19, 4.2e-19, 4e-19,
        3.2e-19, 2.35e-19, 2.15e-19, 1.75e-19, 1.6e-19, 1.54e-19, 1.46e-19,
        1.42e-19, 1.38e-19, 1.33e-19, 1.29e-19, 1.23e-19, 1.18e-19, 1.13e-19,
        1.06e-19, 1.04e-19, 1.02e-19, 1.01e-19, 9.8e-20, 9.1e-20, 8.6e-20,
        7.1e-20, 6.1e-20, 5.5e-20, 4.9e-20, 4.2e-20, 3.3e-20, 2.44e-20,
        1.96e-20, 1.55e-20},
        // H2O
        {1.65e-17, 8.8e-18, 7.7e-18, 7.04e-18, 5.555e-18, 4.73e-18, 3.85e-18,
        3.08e-18, 2.2e-18, 1.43e-18, 7.26e-19, 5.5e-19, 4.4e-19, 3.52e-19,
        2.75e-19, 2.31e-19, 9.9e-20, 6.38e-20, 5.17e-20, 5.28e-20, 5.5e-20,
        5.83e-20, 6.16e-20, 6.435e-20, 6.6e-20, 6.82e-20, 6.875e-20, 6.38e-20,
        6.05e-20, 5.17e-20, 3.74e-20, 2.31e-20, 1.375e-20, 9.9e-21, 6.16e-21},
        // H2
        {6.4e-20, 6.4e-20, 6.5e-20, 6.6e-20, 6.8e-20, 7.1e-20, 7.2e-20,
        7.3e-20, 7.7e-20, 8e-20, 8.5e-20, 8.96e-20, 9.28e-20, 9.85e-20,
        1.05e-19, 1.085e-19, 1.14e-19, 1.16e-19, 1.2e-19, 1.25e-19, 1.3e-19,
        1.345e-19, 1.39e-19, 1.47e-19, 1.63e-19, 1.74e-19, 1.78e-19, 1.8e-19,
        1.825e-19, 1.825e-19, 1.81e-19, 1.79e-19, 1.77e-19, 1.7e-19, 1.64e-19,
        1.6e-19, 1.56e-19, 1.48e-19, 1.4e-19, 1.31e-19, 1.22e-19, 1.04e-19,
        8.9e-20, 7.85e-20, 6e-20, 5.2e-20, 4.5e-20, 4.2e-20, 3.9e-20, 3.6e-20,
        3.4e-20, 2.9e-20, 2.6e-20, 2.3e-20, 1.9e-20, 1.62e-20, 1.28e-20,
        9.2e-21, 7.2e-21, 5.4e-21, 3.7e-21, 2.9e-21, 2.1e-21, 1.4e-21,
        1.04e-21, 7.7e-22},
        // CH4
        {2e-19, 2e-19, 1.16e-19, 7e-20, 4e-20, 2e-20, 1e-20, 6e-21, 7e-21,
        1e-20, 2e-20, 4e-20, 7e-20, 1.18e-19, 1.4e-19, 1.17e-19, 1.12e-19,
        8e-20, 4e-20, 1.18e-20, 1.12e-20, 4e-21, 1e-22}
    };
}

void IonGasTransport::init(thermo_t* thermo, int mode, int log_level)
{
    m_thermo = thermo;
    m_nsp = m_thermo->nSpecies();
    m_mode = mode;
    if (m_mode == CK_Mode) {
        throw CanteraError("IonGasTransport::init(thermo, mode, log_level)",
                           "mode = CK_Mode, which is an outdated lower-order fit.");
    }
    m_log_level = log_level;
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++) {
        if (m_speciesCharge[k] != 0){
            if (k != m_kElectron) {
                m_kIon.push_back(k);
            }
        } else {
            m_kNeutral.push_back(k);
        }
    }
    // set up O2/O2- collision integral [A^2]
    // Data taken from Prager (2005)
    const vector_fp temp{300.0, 400.0, 500.0, 600.0, 800.0, 1000.0,
                         1200.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0};
    const vector_fp om11_O2{120.0, 107.0, 98.1, 92.1, 83.0, 77.0,
                            72.6, 67.9, 62.7, 59.3, 56.7, 53.8};
    vector_fp w(temp.size(),-1);
    int degree = 5;
    m_om11_O2.resize(degree + 1);
    polyfit(temp.size(), degree, temp.data(), om11_O2.data(),
            w.data(), m_om11_O2.data());
    // set up Monchick and Mason parameters
    setupCollisionParameters();
    // set up n64 parameters
    setupN64();
    // setup  collision integrals
    setupCollisionIntegral();
    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);
    m_visc.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_phi.resize(m_nsp, m_nsp, 0.0);
    m_bdiff.resize(m_nsp, m_nsp);
    m_cond.resize(m_nsp);

    // make a local copy of the molecular weights
    m_mw = m_thermo->molecularWeights();

    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t k = j; k < m_nsp; k++) {
            m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
            m_wratjk(k,j) = sqrt(m_wratjk(j,k));
            m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
        }
    }
}

double IonGasTransport::viscosity()
{
    update_T();
    update_C();

    if (m_visc_ok) {
        return m_viscmix;
    }

    double vismix = 0.0;
    // update m_visc and m_phi if necessary
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    multiply(m_phi, m_molefracs.data(), m_spwork.data());

    for (size_t k : m_kNeutral) {
        vismix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    m_viscmix = vismix;
    return vismix;
}

double IonGasTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_spcond_ok) {
        updateCond_T();
    }
    if (!m_condmix_ok) {
        doublereal sum1 = 0.0, sum2 = 0.0;
        for (size_t k : m_kNeutral) {
            sum1 += m_molefracs[k] * m_cond[k];
            sum2 += m_molefracs[k] / m_cond[k];
        }
        m_lambda = 0.5*(sum1 + 1.0/sum2);
        m_condmix_ok = true;
    }
    return m_lambda;
}

double IonGasTransport::electricalConductivity()
{
    vector_fp mobi(m_nsp);
    getMobilities(&mobi[0]);
    double p = m_thermo->pressure();
    double sum = 0.0;
    for (size_t k : m_kIon) {
        double ND_k = m_molefracs[k] * p / m_kbt;
        sum += ND_k * abs(m_speciesCharge[k]) * ElectronCharge * mobi[k];
    }
    if (m_kElectron != npos) {
        sum += m_molefracs[m_kElectron] * p / m_kbt *
               ElectronCharge * mobi[m_kElectron];
    }
    return sum;
}

void IonGasTransport::fitDiffCoeffs(MMCollisionInt& integrals)
{
    GasTransport::fitDiffCoeffs(integrals);

    // number of points to use in generating fit data
    const size_t np = 50;
    int degree = 4;
    double dt = (m_thermo->maxTemp() - m_thermo->minTemp())/(np-1);
    vector_fp tlog(np);
    vector_fp w(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = m_thermo->minTemp() + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1);
    double err = 0.0, relerr = 0.0,
           mxerr = 0.0, mxrelerr = 0.0;

    vector_fp diff(np + 1);
    // The array order still not ideal
    for (size_t k = 0; k < m_nsp; k++) {
        for (size_t j = k; j < m_nsp; j++) {
            if (m_alpha[k] == 0.0 || m_alpha[j] == 0.0 ||
                k == m_kElectron || j == m_kElectron) {
                continue;
            }
            if (m_speciesCharge[k] == 0) {
                if (m_speciesCharge[j] == 0) {
                    continue;
                }
            } else {
                if (m_speciesCharge[j] != 0) {
                    continue;
                }
            }
            for (size_t n = 0; n < np; n++) {
                double t = m_thermo->minTemp() + dt*n;
                double eps = m_epsilon(j,k);
                double tstar = Boltzmann * t/eps;
                double sigma = m_diam(j,k);
                double om11 = omega11_n64(tstar, m_gamma(j,k))
                              * Pi * sigma * sigma;
                // Stockmayer-(n,6,4) model is not suitable for collision
                // between O2/O2- due to resonant charge transfer.
                // Therefore, the experimental collision data is used instead.
                if ((k == m_thermo->speciesIndex("O2-") ||
                     j == m_thermo->speciesIndex("O2-")) &&
                    (k == m_thermo->speciesIndex("O2") ||
                     j == m_thermo->speciesIndex("O2"))) {
                       om11 = poly5(t, m_om11_O2.data()) / 1e20;
                }
                double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/m_reducedMass(k,j))
                    * pow(Boltzmann * t, 1.5) / om11;

                diff[n] = diffcoeff/pow(t, 1.5);
                w[n] = 1.0/(diff[n]*diff[n]);
            }
            polyfit(np, degree, tlog.data(), diff.data(), w.data(), c.data());

            for (size_t n = 0; n < np; n++) {
                double val, fit;
                double t = exp(tlog[n]);
                double pre = pow(t, 1.5);
                val = pre * diff[n];
                fit = pre * poly4(tlog[n], c.data());
                err = fit - val;
                relerr = err/val;
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }
            size_t sum = k * (k + 1) / 2;
            m_diffcoeffs[k*m_nsp+j-sum] = c;
            if (m_log_level >= 2) {
                writelog(m_thermo->speciesName(k) + "__" +
                         m_thermo->speciesName(j) + ": [" + vec2str(c) + "]\n");
            }
        }
    }

    if (m_log_level) {
        writelogf("Maximum binary diffusion coefficient absolute error:"
                 "  %12.6g\n", mxerr);
        writelogf("Maximum binary diffusion coefficient relative error:"
                 "%12.6g", mxrelerr);
    }
}

void IonGasTransport::setupN64()
{
    m_gamma.resize(m_nsp, m_nsp, 0.0);
    for (size_t i : m_kIon) {
        for (size_t j : m_kNeutral) {
            if (m_alpha[j] != 0.0 && m_alpha[i] != 0.0) {
                double r_alpha = m_alpha[i] / m_alpha[j];
                // save a copy of polarizability in Angstrom
                double alphaA_i = m_alpha[i] * 1e30;
                double alphaA_j = m_alpha[j] * 1e30;
                // The ratio of dispersion to induction forces
                double xi = alphaA_i / (m_speciesCharge[i] * m_speciesCharge[i] *
                            (1.0 + pow(2 * r_alpha, 2./3.)) * sqrt(alphaA_j));

                // the collision diameter
                double K1 = 1.767;
                double kappa = 0.095;
                m_diam(i,j) = K1 * (pow(m_alpha[i], 1./3.) + pow(m_alpha[j], 1./3.)) /
                              pow(alphaA_i * alphaA_j * (1.0 + 1.0 / xi), kappa);

                // The original K2 is 0.72, but Han et al. suggested that K2 = 1.44
                // for better fit.
                double K2 = 1.44;
                double epsilon = K2 * ElectronCharge * ElectronCharge *
                                 m_speciesCharge[i] * m_speciesCharge[i] *
                                 m_alpha[j] * (1.0 + xi) /
                                 (8 * Pi * epsilon_0 * pow(m_diam(i,j),4));
                if (epsilon != 0.0) {
                    m_epsilon(i,j) = epsilon;
                }

                // Calculate dipersion coefficient and quadrupole polarizability
                // from curve fitting if not available.
                // Neutrals
                if (m_disp[j] == 0.0) {
                    m_disp[j] = exp(1.8846*log(alphaA_j)-0.4737)* 1e-50;
                }
                if (m_quad_polar[j] == 0.0) {
                    m_quad_polar[j] = 2.0 * m_disp[j];
                }
                // Ions
                if (m_disp[i] == 0.0) {
                    if (m_speciesCharge[i] > 0) {
                        m_disp[i] = exp(1.8853*log(alphaA_i)+0.2682)* 1e-50;
                    } else {
                        m_disp[i] = exp(3.2246*log(alphaA_i)-3.2397)* 1e-50;
                    }
                }

                // The binary dispersion coefficient is determined by the combination rule
                // Reference:
                // Tang, K. T. "Dynamic polarizabilities and van der Waals coefficients."
                // Physical Review 177.1 (1969): 108.
                double C6 = 2.0 * m_disp[i] * m_disp[j] /
                            (1.0/r_alpha * m_disp[i] + r_alpha * m_disp[j]);

                m_gamma(i,j) = (2.0 / pow(m_speciesCharge[i],2) * C6 + m_quad_polar[j]) /
                               (m_alpha[j] * m_diam(i,j) * m_diam(i,j));//Dimensionless

                // properties are symmetric
                m_diam(j,i) = m_diam(i,j);
                m_epsilon(j,i) = m_epsilon(i,j);
                m_gamma(j,i) = m_gamma(i,j);
            }
        }
    }
}

double IonGasTransport::omega11_n64(const double tstar, const double gamma)
{
    double logtstar = log(tstar);
    double om11 = 0.0;
    if (tstar < 0.01) {
        throw CanteraError("IonGasTransport::omega11_n64(tstar, gamma)",
                           "tstar = {} is smaller than 0.01", tstar);
    } else if (tstar <= 0.04) {
        // for interval 0.01 to 0.04, SSE = 0.006; R^2 = 1; RMSE = 0.020
       om11 = 2.97 - 12.0 * gamma
              - 0.887 * logtstar
              + 3.86 * gamma * gamma
              - 6.45 * gamma * logtstar
              - 0.275 * logtstar * logtstar
              + 1.20 * gamma * gamma * logtstar
              - 1.24 * gamma * logtstar * logtstar
              - 0.164 * pow(logtstar,3);
    } else if (tstar <= 1000) {
        // for interval 0.04 to 1000, SSE = 0.282; R^2 = 1; RMSE = 0.033
       om11 = 1.22 - 0.0343 * gamma
              + (-0.769 + 0.232 * gamma) * logtstar
              + (0.306 - 0.165 * gamma) * logtstar * logtstar
              + (-0.0465 + 0.0388 * gamma) * pow(logtstar,3)
              + (0.000614 - 0.00285 * gamma) * pow(logtstar,4)
              + 0.000238 * pow(logtstar,5);
    } else {
        throw CanteraError("IonGasTransport::omega11_n64(tstar, gamma)",
                           "tstar = {} is larger than 1000", tstar);
    }
    return om11;
}

void IonGasTransport::getMixDiffCoeffs(double* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    double mmw = m_thermo->meanMolecularWeight();
    double p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            if (k == m_kElectron) {
                d[k] = 0.4 * m_kbt / ElectronCharge;
            } else {
                double sum2 = 0.0;
                for (size_t j : m_kNeutral) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
                if (sum2 <= 0.0) {
                    d[k] = m_bdiff(k,k) / p;
                } else {
                    d[k] = (mmw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
                }
            }
        }
    }
}

void IonGasTransport::getMobilities(double* const mobi)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    double p = m_thermo->pressure();
    for (size_t k = 0; k < m_nsp; k++) {
        if (k == m_kElectron) {
            mobi[k] = 0.4;
        } else {
            mobi[k] = 0.0;
        }
    }
    for (size_t k : m_kIon) {
        double sum = 0.0;
        for (size_t j : m_kNeutral) {
            double bmobi = m_bdiff(k,j) * ElectronCharge / m_kbt;
            sum += m_molefracs[j] / bmobi;
        }
        mobi[k] = 1.0 / sum / p;
    }
}

}
