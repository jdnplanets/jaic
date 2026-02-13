/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // processElec.h

#pragma once

#include "device_arrays.h"
#include "Precip.hpp"


// Function declarations
void runPrimariesKernel(DeviceArrays darrs, DeviceArrays harrs, SimParams p);
