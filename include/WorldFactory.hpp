/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // WorldFactory.hpp

#pragma once
#include <memory>
#include <string>
#include "WorldOverrides.hpp"

struct World1D;

std::unique_ptr<World1D> make_world(const std::string& name,
                                   const WorldOverrides& overrides = {});
