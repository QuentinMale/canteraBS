//! @file findinputfile.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_FINDINPUTFILE_H
#define CT_FINDINPUTFILE_H

#include "application.h"

extern "C" {
    void findInputFileChar(const char* name, const char* fileName) {
        const std::string nameString = std::string(name);
        fileName = Cantera::Application::findInputFile(nameString).c_str();
    };
}
