// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/Collision.h"

namespace Cantera {

Collision::Collision()
    : threshold(0.0)
{
}

void Collision::setParameters(const AnyMap& node)
{
    if (node.empty()) {
        return;
    }
    equation = node["equation"].asString();
    kind = node["kind"].asString();
    threshold = node.getDouble("threshold", 0.0);
    AnyValue dataNode = node["data"];
    energy_data = dataNode["energy"].asVector<double>();
    cross_section_data = dataNode["cross-section"].asVector<double>();
    id = node["id"].asString();
    input = node;
}

void Collision::getParameters(AnyMap& collisionNode) const
{
    collisionNode["equation"] = equation;
    collisionNode["kind"] = kind;
    collisionNode["threshold"] = threshold;
    collisionNode["energy-data"] = energy_data;
    collisionNode["cross-section-data"] = cross_section_data;
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

unique_ptr<Collision> newCollision(const AnyMap& node, const Kinetics& kin)
{
    unique_ptr<Collision> ec(new Collision());
    ec->setParameters(node);
    return ec;
}

std::vector<shared_ptr<Collision>> getCollisions(const AnyValue& items, const Kinetics& kin)
{
    std::vector<shared_ptr<Collision> > all_collisions;
    for (const auto& node : items.asVector<AnyMap>()) {
        all_collisions.emplace_back(newCollision(node, kin));
    }
    return all_collisions;
}

}
