/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
 This file implements the make_world factory function, which constructs and initializes World1D 
 */

 // WorldFactory.cpp

#include "WorldFactory.hpp"
#include "World1D.hpp"
#include "BrownDwarf.hpp"
#include "Jupiter.hpp"

#include "NullIonChemistryModel.hpp"
#include "H3CH5IonChemistryModel.hpp"

#include <stdexcept>
#include <cstdlib>   // getenv

static std::string home_dir() {
    const char* h = std::getenv("HOME");
    return h ? std::string(h) : std::string();
}

std::unique_ptr<World1D> make_world(const std::string& name,
                                   const WorldOverrides& overrides)
{
    std::unique_ptr<World1D> w;

    if (name == "BrownDwarf" || name == "Brown Dwarf") {
        w = std::make_unique<BrownDwarf>(BrownDwarf::Defaults());

        // No chemistry: all ion production becomes H3+.
        w->ionChemistryModel = std::make_shared<NullIonChemistryModel>();

    } else if (name == "Jupiter") {
        w = std::make_unique<Jupiter>(); // defaults inside class

        // H3+ <-> CH5+ chemistry (steady-state, no transport)
        w->ionChemistryModel = std::make_shared<H3CH5IonChemistryModel>();

    } else {
        throw std::invalid_argument("Unknown world: " + name);
    }

    w->apply_overrides(overrides); // set params before init()
    w->init();                     // populate arrays consistently
    return w;
}
