/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // SpeciesData.hpp

#pragma once
#include "constants.h"
#include <vector>

using namespace constants;

struct Species {
    std::string name;
    double alpha;
    double r;
    double m;   // in kg
    double q; // in C
};

inline Species Electron{
    "e-",              // name
    1,                  // alpha in m^2
    2.82e-15,          // radius in m (classical electron radius)
    m_e,               // mass in kg
    -qe                // charge in C
};

inline Species H3p{
    "H3+",             // name
    1.15e-7,            // alpha in m^2
    2*0.52e-10,           // radius in m(bond length/sqrt(3) from geometry)
    3 * m_p,            // mass in kg
    qe                  // charge in C
};

inline Species CH5p{
    "CH5+",            // name
    2.7e-7,            // alpha in m^2
    1.1e-10,           // radius in m(bond length/sqrt(3) from geometry)
    17 * m_p,            // mass in kg
    qe                  // charge in C
};

inline Species C3Hnp{
    "C3Hn+",           // name
    7.5e-7,            // alpha in m^2
    3 * 1.1e-10,           // radius in m(bond length/sqrt(3) from geometry)
    40 * m_p,            // mass in kg
    qe                  // charge in C
};


inline Species H2{
    "H2",              // name
    -1,                  // alpha in m^2
    2*0.74e-10,          // H2 radius in m (bond length)
    2 * m_p,            // mass in kg
    0                   // neutral species has no charge
};