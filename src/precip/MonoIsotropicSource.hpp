/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
 This file defines a MonoIsotropicSource class, which is a specific implementation of the abstract Source class. 
 The MonoIsotropicSource generates a monenergetic, isotropic distribution of particles for the precipitation simulation. 
 */

 // MonoIsotropicSource.hpp

#pragma once
#include "Source.hpp"
#include "constants.h"
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace constants;
using namespace std;
#include "Rnd.h"


class MonoIsotropicSource : public Source {
public:
  void init(
    const Params& sp,
    std::vector<float>& z,
    std::vector<float>& y,
    std::vector<float>& vz,
    std::vector<float>& vy,
    std::vector<float>& E,
    std::vector<int>&   alive
  ) const override {
    vz.resize(sp.N);
    vy.resize(sp.N);
    float v0 = Etov(sp.Einit); // Convert energy in eV to velocity in m/s
    for (size_t i = 0; i < sp.N; ++i) {
        float theta = pi * rnd() - pi/2; // uniform distribution
        vz[i] = -v0 * std::cos(theta);
        vy[i] = v0 * std::sin(theta);
    }
    z.assign(sp.N, sp.Zinit);
    y.assign(sp.N, 0.0f);
    E.assign(sp.N, sp.Einit); // in eV
    alive.assign(sp.N, 1);
  
  }

    /// return a filename‐friendly label for this source
    string label(const Params& sp) const override {
    ostringstream os;
    if (sp.Einit < 1e3) {
      os << "mono_" << static_cast<int>(sp.Einit) << "eV";
    } else  
    os << "mono_" << setfill('0') << setw(3) << static_cast<int>(sp.Einit/1e3) << "keV";
    return os.str();
  }
  
  string banner(const Params& sp) const override {
    ostringstream os;
    if (sp.Einit < 1e3) {
      os << "Monenergetic Isotropic beam @ " << sp.Einit << " eV\n";
    } else  
    os << "Monenergetic Isotropic beam @ " << sp.Einit/1e3 << " keV\n";
    return os.str();
  }
};

