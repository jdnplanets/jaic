/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // constants.h

//constants.hpp
#pragma once 
#ifdef USE_CUDA
#include <cuda_runtime.h> // for CUDA runtime API
#endif

namespace constants {
#ifdef USE_CUDA
    // floats because of CUDA compatibility
    __constant__ constexpr float pi = 3.14159265358979323846;
    __constant__ constexpr float m_e = 9.10938356e-31;    // kg - mass of an electron
    __constant__ constexpr float m_p = 1.6726219e-27;    // kg - mass of a proton
    __constant__ constexpr float eV = 1.60217662e-19;    // J - electron volt
    __constant__ constexpr float qe = 1.60217662e-19;     // C - elementary charge
    __constant__ constexpr float e = 1.60217662e-19;     // C - elementary charge
    __constant__ constexpr float c = 299792458;          // m/s - speed of light
    __constant__ constexpr float k_B = 1.380649e-23;     // J/K - Boltzmann constant
    __constant__ constexpr float h = 6.62607015e-34;    // J s - Planck constant
#else
    constexpr float pi = 3.14159265358979323846f;
    constexpr float m_e = 9.10938356e-31f;    // kg - mass of an electron
    constexpr float m_p = 1.6726219e-27f;    // kg - mass of a proton
    constexpr float eV = 1.60217662e-19f;    // J - electron volt
    constexpr float qe = 1.60217662e-19f;     // C - elementary charge
    constexpr float e = 1.60217662e-19f;     // C - elementary charge
    constexpr float c = 299792458.0f;        // m/s - speed of light
    constexpr float k_B = 1.380649e-23f;     // J/K - Boltzmann constant
    constexpr float h = 6.62607015e-34f;    // J s - Planck constant
#endif

}