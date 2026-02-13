/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // IonChemistryRates.hpp

#pragma once

// Central location for chemistry rate coefficients.

#include <algorithm>

#include "SpeciesData.hpp"

namespace IonChemistryRates {

    // Reaction: H3+ + CH4 -> CH5+ + ...
    // Units: cm^3 s^-1
    inline constexpr double k_H3p_CH4 = 2.0e-9;

    // Dissociative recombination coefficients.
    // Units: cm^3 s^-1
    //
    inline double alpha_H3p(double T_K) {
        (void)T_K;
        return 1.15e-7*std::pow(300./T_K, 0.65);
    }

    inline double alpha_CH5p(double T_K) {
        (void)T_K;
        return 2.78e-7*std::pow(300./T_K, 0.52);
    }

    inline double alpha_C3Hnp(double T_K) {
        (void)T_K;
        return 3.0e-7*std::pow(300./T_K, 0.50);
    }

}
