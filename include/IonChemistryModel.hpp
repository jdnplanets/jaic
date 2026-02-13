/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the IIonChemistryModel interface and a specific implementation 
(H3CH5IonChemistryModel) for a simple two-ion chemistry system. 
The IIonChemistryModel interface allows for different ion chemistry models to be used in the World1D class, 
which represents the 1D atmospheric model. 

*/

 // IonChemistryModel.hpp

#pragma once
#include <vector>
#include "World1D.hpp"

// Forward-declare Params
struct Params;

// Interface for ion chemistry models (steady-state, no transport).
// A World chooses an implementation in WorldFactory.
struct IIonChemistryModel {
    virtual ~IIonChemistryModel() = default;

    // Qtot[z] is the total ion production rate at altitude z (cm^-3 s^-1).
    // Implementations write results into world.nH3p/nCH5p/nC3Hnp/world.ne.
    virtual void compute(const Params& sp,
                         World1D& world,
                         const std::vector<double>& Qtot) const = 0;
};
