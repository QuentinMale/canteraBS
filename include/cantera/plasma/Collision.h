//! @file Collision.h Declaration for class Cantera::Collision.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_Collision_H
#define CT_Collision_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
class ThermoPhase;

//! Contains data about the cross sections of collision
/*!
 *  This class stores the cross section data of collision
 */
class Collision
{
public:
    Collision();

    //! Collision objects are not copyable or assignable
    Collision(const Collision&) = delete;
    Collision& operator=(const Collision& other) = delete;
    ~Collision();

    //! Return the parameters such that an identical Collision could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    //! @param withInput  If true, include additional input data fields associated
    //!   with the object, such as user-defined fields from a YAML input file, as
    //!   contained in the #input attribute.
    AnyMap parameters(bool withInput=true) const;

    //! Set up reaction based on AnyMap *node*
    virtual void setParameters(const AnyMap& node);

    //! The name of the kind of collision
    std::string kind;

    //! True if the current reaction is marked as duplicate
    bool duplicate;

    //! True if the current reaction is reversible. False otherwise
    bool reversible;

    //! Reactant species and stoichiometric coefficients
    Composition reactants;

    //! Product species and stoichiometric coefficients
    Composition products;

    //! Energy data
    vector_fp energy_data;

    //! Cross section data
    vector_fp cross_section_data;

    //! The threshold of a process
    double threshold;

    //! The chemical equation for this reaction
    std::string equation() const;

    //! The reactant side of the chemical equation for this reaction
    virtual std::string reactantString() const;

    //! The product side of the chemical equation for this reaction
    virtual std::string productString() const;

    //! Input parameters used to define a collision, e.g. from a YAML input file.
    AnyMap input;

protected:
    //! Store the parameters of a Reaction needed to reconstruct an identical
    //! object using the newReaction(AnyMap&, Kinetics&) function. Does not
    //! include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& collisionNode) const;
};

//! Parse collision equation
void parseCollisionEquation(Collision& C, const AnyValue& equation);

//! create an Collision object to store data.
unique_ptr<Collision> newCollision(const AnyMap& node);

//! Get a vector of Collision objects to access the data.
std::vector<shared_ptr<Collision>> getCollisions(const AnyValue& items);

}

#endif
