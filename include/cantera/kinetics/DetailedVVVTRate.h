//! @file DetailedVVVTRate.h   Header for plasma reaction rates corresponding to detailed vibration handling).

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DETAILEDVVVTRATE_H
#define CT_DETAILEDVVVTRATE_H

#include "Arrhenius.h"

namespace Cantera
{

//! Data container holding shared data specific to TwoTempPlasmaRate
/**
 * The data container `TwoTempPlasmaData` holds precalculated data common to
 * all `TwoTempPlasmaRate` objects.
 */
struct DetailedVibData : public ReactionData
{
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
};


class DetailedVVVTRate : public ArrheniusBase
{
public:
    DetailedVVVTRate();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param scaling  The scaling factor from the harmonic oscillator theory
     *  @param B  The constant in the exponential
     *  @param C  The constant over T**1/3 in the exponential
     *  @param D  The constant over T**2/3 in the exponential
     *  @param alpha The exponent of the pre-exponential temperature
     */
    DetailedVVVTRate(double A, double B, double C, double D, double b, double scaling);

    explicit DetailedVVVTRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<DetailedVVVTRate, DetailedVibData>>();
    }

    const string type() const override {
        return "detailed-vv-vt";
    }

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const DetailedVibData& shared_data) const {
        return m_scaling * m_A * std::exp(m_b * shared_data.logT + m_B + m_C * std::pow(shared_data.recipT, 1.0/3.0) + m_D * std::pow(shared_data.recipT, 2.0/3.0));
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  This method does not consider changes of electron temperature.
     *  A corresponding warning is raised.
     *  @param shared_data  data shared by all reactions of a given type
     */
    double ddTScaledFromStruct(const DetailedVibData& shared_data) const;


};

}

#endif
