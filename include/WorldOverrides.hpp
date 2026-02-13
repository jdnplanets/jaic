/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // WorldOverrides.hpp

#pragma once
#include <optional>

// A generic “override bag”.
struct WorldOverrides {
    
    std::optional<double> B_T;
    std::optional<double> Tiso_K;   // for isothermal worlds
    std::optional<double> g_m_s2;   // gravity
    std::optional<double> n0_cm3;   // reference density
    std::optional<double> Zref_m;   // reference altitude for n0

};
