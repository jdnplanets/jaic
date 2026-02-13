/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the MonoSource class, which is a specific implementation of the Source interface for 
generating a monoenergetic beam of particles. The MonoSource class initializes the particle arrays with a 
specified energy and vertical velocity.
*/

 // MonoSource.hpp

#pragma once
#include "Source.hpp"
#include "constants.h"
#include <sstream>
#include <cmath>
#include <iomanip> // For std::setprecision
// using namespace constants;
using namespace std;


class MonoSource : public Source {
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
    float v0 = Etov(sp.Einit); // Convert energy in eV to velocity in m/s
    z.assign(sp.N, sp.Zinit);
    y.assign(sp.N, 0.0f);
    vz.assign(sp.N, -v0);
    vy.assign(sp.N, 0.0f);
    E.assign(sp.N, sp.Einit); // in eV
    alive.assign(sp.N, 1);
  }

    /// return a filename‐friendly label for this source
    string label(const Params& sp) const override {
    ostringstream os;
    if (sp.Einit < 1e3) {
      int int_part = static_cast<int>(sp.Einit);
      int frac_part = static_cast<int>(round((sp.Einit - int_part) * 100));
      os << "mono_" << setfill('0') << setw(4) << int_part << "_"
         << setw(2) << frac_part << "eV";
        } else {
      double val = sp.Einit / 1e3;
      int int_part = static_cast<int>(val);
      int frac_part = static_cast<int>(round((val - int_part) * 100));
      os << "mono_" << setfill('0') << setw(4) << int_part << "_" 
        << setw(2) << frac_part << "keV";
        }
    return os.str();
  }
    
  string banner(const Params& sp) const override {
    ostringstream os;
    if (sp.Einit < 1e3) {
      os << "Monenergetic vertical beam @ " << sp.Einit << " eV\n";
    } else  
    os << "Monenergetic vertical beam @ " << sp.Einit/1e3 << " keV\n";
    return os.str();
  }
};

