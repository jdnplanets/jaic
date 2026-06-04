/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the DeviceArrays struct, which contains pointers to arrays used in the GPU computations 
for the precipitation simulation. These arrays include lookup tables for cross sections and probabilities, 
as well as state variables for the particles being simulated.
*/

 // device_arrays.h

#pragma once


struct DeviceArrays {

    // Lookup tables
    float* Zinx;
    float* nH2; // neutral H2 density in cm^-3
    float* sigmasE;
    float* total_sigma;
    float* sigmas;
    float* prob;
    int* alias;
    float* deltaEs;

    //variables
    float* dt;
    int* zcell;
    float* zlocal;
    float* y;
    float* vz;
    float* vy;
    float* E; // energy in eV
    int* alive;
    int* ps; //  1 for primary, 2 for secondary
    // float* zfinal;
    // float* yfinal;
    int* nion;
    int* nion2;
    int* theta_sampled;
    int* colcount;
    int* nexB;
    int* nexC;
    int* nexEF; // number of excitations to singlet E

};