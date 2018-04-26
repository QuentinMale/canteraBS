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
}

void IonGasTransport::setupMM()
{
    GasTransport::setupMM();
    m_gamma.resize(m_nsp, m_nsp, 0.0);
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            if (k != m_kElectron) {
                m_kCharge.push_back(k);
            }
        } else {
            m_kNeutral.push_back(k);
        }
    }
    setupN64();
}

void IonGasTransport::fitDiffCoeffs(MMCollisionInt& integrals)
{
    GasTransport::fitDiffCoeffs(integrals);

    // number of points to use in generating fit data
    const size_t np = 50;
    int degree = (m_mode == CK_Mode ? 3 : 4);
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
    for (size_t k : m_kNeutral) {
        for (size_t j : m_kCharge) {
            if (j >= k && m_alpha[k] != 0.0) {
                size_t sum = 0;
                for (size_t i = 0; i <= k; i++) {
                    sum += i;
                }
                for (size_t n = 0; n < np; n++) {
                    double t = m_thermo->minTemp() + dt*n;
                    double eps = m_epsilon(j,k);
                    double tstar = Boltzmann * t/eps;
                    double sigma = m_diam(j,k);
                    double om11 = omega11_n64(tstar, m_gamma(j,k));
                    double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/m_reducedMass(k,j))
                        * pow(Boltzmann * t, 1.5) / (Pi * sigma * sigma * om11);

                    if (m_mode == CK_Mode) {
                        diff[n] = log(diffcoeff);
                        w[n] = -1.0;
                    } else {
                        diff[n] = diffcoeff/pow(t, 1.5);
                        w[n] = 1.0/(diff[n]*diff[n]);
                    }
                }
                polyfit(np, degree, tlog.data(), diff.data(), w.data(), c.data());

                for (size_t n = 0; n < np; n++) {
                    double val, fit;
                    if (m_mode == CK_Mode) {
                        val = exp(diff[n]);
                        fit = exp(poly3(tlog[n], c.data()));
                    } else {
                        double t = exp(tlog[n]);
                        double pre = pow(t, 1.5);
                        val = pre * diff[n];
                        fit = pre * poly4(tlog[n], c.data());
                    }
                    err = fit - val;
                    relerr = err/val;
                    mxerr = std::max(mxerr, fabs(err));
                    mxrelerr = std::max(mxrelerr, fabs(relerr));
                }
                m_diffcoeffs[k*m_nsp+j-sum] = c;
                if (m_log_level >= 2) {
                    writelog(m_thermo->speciesName(k) + "__" +
                             m_thermo->speciesName(j) + ": [" + vec2str(c) + "]\n");
                }
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
    // Reference:
    // Selle, Stefan, and Uwe Riedel.
    // "Transport coefficients of reacting air at high temperatures."
    // AIAA 211 (2000): 10-13.
    // The method of estimating the potential parameters is from:
    // Aquilanti, Vincenzo, David Cappelletti, and Fernando Pirani.
    // "Range and strength of interatomic forces: dispersion and induction
    // contributions to the bonds of dications and of ionic molecules."
    // Chemical physics 209.2 (1996): 299-311.

    for (size_t i : m_kCharge) {
        for (size_t j : m_kNeutral) {
            if (m_alpha[j] != 0.0) {
                double r_alpha = m_alpha[i] / m_alpha[j];
                // save a copy of polarizability in Angstrom
                double alphaA_i = m_alpha[i] * 1e30;
                double alphaA_j = m_alpha[j] * 1e30;
                // The ratio of dispersion to induction forces
                double xi = m_speciesCharge[i] * m_speciesCharge[i];
                       xi *= 1.0 + pow((2 * r_alpha ),(2./3.));
                       xi *= sqrt(alphaA_j);
                       xi = alphaA_i / xi;

                // the collision diameter
                m_diam(i,j) = 1.767;
                m_diam(i,j) *= pow(m_alpha[i],(1./3.)) + pow(m_alpha[j],(1./3.));
                m_diam(i,j) /= pow((alphaA_i * alphaA_j * (1.0 + 1.0 / xi)),0.0095);

                double epsilon = 0.72 * ElectronCharge * ElectronCharge;
                epsilon *= m_speciesCharge[i] * m_speciesCharge[i];
                epsilon *= m_alpha[j] * (1.0 + xi);
                epsilon /= 8 * Pi * epsilon_0 * pow(m_diam(i,j),4);

                if (epsilon != 0.0) {
                    m_epsilon(i,j) = epsilon;
                }

                // Calculate dipersion coefficient and quadrupole polarizability
                // from curve fitting if not available.
                // Reference:
                // Han, Jie, et al. "Numerical modelling of ion transport in flames."
                // Combustion Theory and Modelling 19.6 (2015): 744-772.
                // Neutrals
                if (m_disp[j] == 0.0) {
                    m_disp[j] = exp(1.8846*log(alphaA_j)-0.4737)* 1e-50;
                }
                if (m_qua_polar[j] == 0.0) {
                    m_qua_polar[j] = 2.0 * m_disp[j];
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

                m_gamma(i,j) = 2.0 / pow(m_speciesCharge[i],2) * C6 + m_qua_polar[j];
                m_gamma(i,j) /= m_alpha[j] * m_diam(i,j) * m_diam(i,j);//Dimensionless

                // properties are symmetric (not sure)
                m_diam(j,i) = m_diam(i,j);
                m_epsilon(j,i) = m_epsilon(i,j);
                m_gamma(j,i) = m_gamma(i,j);
            }
        }
    }
}

double IonGasTransport::omega11_n64(const double tstar, const double gamma)
{
    // collision integral should be able to use proper table to do the polyfit
    // Reference:
    // Viehland, L. A., et al. "Tables of transport collision integrals for 
    // (n, 6, 4) ion-neutral potentials." 
    // Atomic Data and Nuclear Data Tables 16.6 (1975): 495-514.
    // However, a known polyfit is used from
    // Han, Jie, et al. "Numerical modelling of ion transport in flames." 
    // Combustion Theory and Modelling 19.6 (2015): 744-772.
    // Note: Han release the range to 1000, but Selle suggested to 10
    double logtstar = log(tstar);
    double om11 = 0.0;
    if (tstar < 0.01) {
        throw CanteraError("IonGasTransport::omega11_n64(tstar, gamma)",
                           "tstar = ", tstar," is smaller than 0.01");
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
                           "tstar = ", tstar, " is larger than 1000");
    }
    return om11;
}

}

