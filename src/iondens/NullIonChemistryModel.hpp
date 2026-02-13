/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
// Simplest chemistry: all ion production becomes H3+; other ions are zero.
*/

 // NullIonChemistryModel.hpp

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "IonChemistryModel.hpp"
#include "IonChemistryRates.hpp"


struct NullIonChemistryModel final : IIonChemistryModel {
    void compute(const Params&,
                 World1D& world,
                 const std::vector<double>& Qtot) const override
    {
        for (int z = 0; z < world.nz; ++z) {
            const double q = std::max(0.0, Qtot[z]); // cm^-3 s^-1
            const double T = (z < static_cast<int>(world.T.size())) ? world.T[z] : 0.0;
            const double alpha = IonChemistryRates::alpha_H3p(T); // cm^3 s^-1

            const double nH3 = (q > 0.0 && alpha > 0.0) ? std::sqrt(q / alpha) : 0.0; // cm^-3

            world.nH3p[z]   = static_cast<float>(nH3);
            world.nCH5p[z]  = 0.0f;
            world.nC3Hnp[z] = 0.0f;
            world.ne[z]     = world.nH3p[z];
        }
    }
};
