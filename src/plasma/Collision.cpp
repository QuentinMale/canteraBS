// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/Collision.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/Array.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"
#include <sstream>
#include <set>

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace Cantera {

Collision::Collision()
    : threshold(0.0)
{
}

Collision::~Collision()
{
}

void Collision::setParameters(const AnyMap& node)
{
    if (node.empty()) {
        return;
    }
    parseCollisionEquation(*this, node["equation"]);
    kind = node["kind"].asString();
    threshold = node.getDouble("threshold", 0.0);
    duplicate = node.getBool("duplicate", false);
    AnyValue dataNode = node["data"];
    energy_data = dataNode["energy"].asVector<double>();
    cross_section_data = dataNode["cross-section"].asVector<double>();
    input = node;
}

void Collision::getParameters(AnyMap& collisionNode) const
{
    collisionNode["equation"] = equation();
    collisionNode["kind"] = kind;
    collisionNode["threshold"] = threshold;
    collisionNode["energy-data"] = energy_data;
    collisionNode["cross-section-data"] = cross_section_data;
    if (duplicate) {
        collisionNode["duplicate"] = true;
    }
}

AnyMap Collision::parameters(bool withInput) const
{
    AnyMap out;
    getParameters(out);
    if (withInput) {
        out.update(input);
    }

    static bool reg = AnyMap::addOrderingRules("Collision",
        {{"head", "kind"},
         {"head", "equation"},
         {"tail", "duplicate"},
        });
    if (reg) {
        out["__type__"] = "Collision";
    }
    return out;
}

std::string Collision::equation() const
{
    if (reversible) {
        return reactantString() + " <=> " + productString();
    } else {
        return reactantString() + " => " + productString();
    }
}

std::string Collision::reactantString() const
{
    std::ostringstream result;
    for (auto iter = reactants.begin(); iter != reactants.end(); ++iter) {
        if (iter != reactants.begin()) {
            result << " + ";
        }
        if (iter->second != 1.0) {
            result << iter->second << " ";
        }
        result << iter->first;
    }
  return result.str();
}

std::string Collision::productString() const
{
    std::ostringstream result;
    for (auto iter = products.begin(); iter != products.end(); ++iter) {
        if (iter != products.begin()) {
            result << " + ";
        }
        if (iter->second != 1.0) {
            result << iter->second << " ";
        }
        result << iter->first;
    }
  return result.str();
}

void parseCollisionEquation(Collision& C, const AnyValue& equation)
{
    // Parse the reaction equation to determine participating species and
    // stoichiometric coefficients
    std::vector<std::string> tokens;
    tokenizeString(equation.asString(), tokens);
    tokens.push_back("+"); // makes parsing last species not a special case

    size_t last_used = npos; // index of last-used token
    bool reactants = true;
    for (size_t i = 1; i < tokens.size(); i++) {
        if (tokens[i] == "+" || ba::starts_with(tokens[i], "(+") ||
            tokens[i] == "<=>" || tokens[i] == "=" || tokens[i] == "=>") {
            std::string species = tokens[i-1];

            double stoich;
            if (last_used != npos && tokens[last_used] == "(+") {
                // Falloff third body with space, e.g. "(+ M)"
                species = "(+" + species;
                stoich = -1;
            } else if (last_used == i-1 && ba::starts_with(species, "(+")
                       && ba::ends_with(species, ")")) {
                // Falloff 3rd body written without space, e.g. "(+M)"
                stoich = -1;
            } else if (last_used == i-2) { // Species with no stoich. coefficient
                stoich = 1.0;
            } else if (last_used == i-3) { // Stoich. coefficient and species
                try {
                    stoich = fpValueCheck(tokens[i-2]);
                } catch (CanteraError& err) {
                    throw InputFileError("parseReactionEquation", equation,
                        err.getMessage());
                }
            } else {
                throw InputFileError("parseReactionEquation", equation,
                    "Error parsing reaction string '{}'.\n"
                    "Current token: '{}'\nlast_used: '{}'",
                    equation.asString(), tokens[i],
                    (last_used == npos) ? "n/a" : tokens[last_used]);
            }

            if (reactants) {
                C.reactants[species] += stoich;
            } else {
                C.products[species] += stoich;
            }

            last_used = i;
        }

        // Tokens after this point are part of the products string
        if (tokens[i] == "<=>" || tokens[i] == "=") {
            C.reversible = true;
            reactants = false;
        } else if (tokens[i] == "=>") {
            C.reversible = false;
            reactants = false;
        }
    }
}

unique_ptr<Collision> newCollision(const AnyMap& node)
{
    unique_ptr<Collision> ec(new Collision());
    ec->setParameters(node);
    return ec;
}

std::vector<shared_ptr<Collision>> getCollisions(const AnyValue& items)
{
    std::vector<shared_ptr<Collision> > all_collisions;
    for (const auto& node : items.asVector<AnyMap>()) {
        all_collisions.emplace_back(newCollision(node));
    }
    return all_collisions;
}

}
