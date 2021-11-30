//! @file Collision.h Declaration for class Cantera::Collision.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_Collision_H
#define CT_Collision_H

#include "cantera/kinetics/Reaction.h"

namespace Cantera
{
class ThermoPhase;
class Kinetics;

//! Contains data about the cross sections of collision
/*!
 *  This class stores the cross section data of collision
 */
class Collision
{
public:
    Collision();

    Collision(const AnyMap& node, const Kinetics& kin);

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

    //! Energy data
    vector_fp energy_data;

    //! Cross section data
    vector_fp cross_section_data;

    //! The threshold of a process
    double threshold;

    //! The chemical equation for this reaction
    std::string equation;

    //! An identification string for the collision, used to link to a reaction
    std::string id;

    //! Input parameters used to define a collision, e.g. from a YAML input file.
    AnyMap input;

protected:
    //! Store the parameters of a Reaction needed to reconstruct an identical
    //! object using the newReaction(AnyMap&, Kinetics&) function. Does not
    //! include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& collisionNode) const;
};

//! create an Collision object to store data.
unique_ptr<Collision> newCollision(const AnyMap& node);

//! Get a vector of Collision objects to access the data.
std::vector<shared_ptr<Collision>> getCollisions(const AnyValue& items);

}

#endif
