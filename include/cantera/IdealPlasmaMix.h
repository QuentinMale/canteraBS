//! @file IdealPlasmaMix.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CXX_IDEALPLASMAMIX
#define CXX_IDEALPLASMAMIX

#include "thermo/IdealGasPhase.h"
#include "kinetics/PlasmaKinetics.h"
#include "kinetics/importKinetics.h"
#include "base/stringUtils.h"

namespace Cantera
{

//! Convenience class which inherits from both IdealGasPhase and GasKinetics
class IdealPlasmaMix :
    public IdealGasPhase,
    public PlasmaKinetics
{
public:
    IdealPlasmaMix() : m_ok(false), m_r(0) {}

    IdealPlasmaMix(const std::string& infile, std::string id_="") :
        m_ok(false), m_r(0)
    {
        m_r = get_XML_File(infile);
        m_id = id_;
        if (id_ == "-") {
            id_ = "";
        }
        m_ok = buildSolutionFromXML(*m_r,
                                    m_id, "phase", this, this);
        if (!m_ok) throw CanteraError("IdealPlasmaMix",
                                          "Cantera::buildSolutionFromXML returned false");
    }

    IdealPlasmaMix(XML_Node& root,
                std::string id_) : m_ok(false), m_r(&root), m_id(id_) {
        m_ok = buildSolutionFromXML(root, id_, "phase", this, this);
    }

    IdealPlasmaMix(const IdealPlasmaMix& other) : m_ok(false),
        m_r(other.m_r),
        m_id(other.m_id) {
        m_ok = buildSolutionFromXML(*m_r, m_id, "phase", this, this);
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }
    friend std::ostream& operator<<(std::ostream& s, IdealPlasmaMix& mix) {
        std::string r = mix.report(true);
        s << r;
        return s;
    }

protected:
    bool m_ok;
    XML_Node* m_r;
    std::string m_id;
};
}

#endif
