/**
 *  @file ElectronCrossSection.cpp Definition file for class ElectronCrossSection.
 */
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCrossSection.h"
#include "cantera/base/global.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : threshold(0.0)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

// void ElectronCrossSection::validate()
// {
//     if (kind == "ionization" || kind == "attachment" || kind == "excitation") {
//         if (threshold < 0.0) {
//             throw CanteraError("ElectronCrossSection::validate",
//                                "The threshold of the process",
//                                "(kind = '{}', target = '{}', product = '{}')",
//                                "cannot be negative", kind, target, product);
//         }
//     } else if (kind != "effective" && kind != "elastic") {
//         throw CanteraError("ElectronCrossSection::validate",
//             "'{}' is an unknown type of cross section data.", kind);
//     }
// }

unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node)
{
    unique_ptr<ElectronCrossSection> ecs(new ElectronCrossSection());

    // XML_Node& dataNode = node.child("data");
    // node["energy-levels"].asVector<double>()

    ecs->kind = node["kind"].asString();
    ecs->target = node["target"].asString();

    vector<double> coeffs_flat;
    //getFloatArray(node, coeffs_flat, false,"", "data");
    coeffs_flat = node["data"].asVector<double>();
    //std::vector<vector_fp> data = node["data"].asVector<vector_fp>();
    // transpose data
    for (size_t i = 0; i < coeffs_flat.size()/2; i++) {
        ecs->energyLevel.push_back(coeffs_flat[2*i]);
        ecs->crossSection.push_back(coeffs_flat[2*i+1]);
    }

    if (node.hasKey("threshold")){
    	ecs->threshold = node["threshold"].asDouble(); //std::stof(node.attrib("threshold"));
    } else {
    	ecs->threshold = 0.0;
    }

    if (node.hasKey("product")) {
    	ecs->product = node["product"].asString();
    } else {
    	ecs->product = ecs->target;
    }

    // Some writelog to check the datas loaded concerning the cross section
    // writelog("Kind :    {:s}\n",ecs->kind);
    // writelog("Target :    {:s}\n",ecs->target);
    // writelog("Product :    {:s}\n",ecs->product);
    // writelog("Threshold :    {:14.5g} eV\n",ecs->threshold);
    // writelog("Energy :  \n");
    // for (size_t i = 0; i < ecs->energyLevel.size(); i++){
    // 	writelog("{:7.2g} {:7.2g} \n",ecs->energyLevel[i], ecs->crossSection[i]);
    // }

    return ecs;
}

}