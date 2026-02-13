/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the UniformSource class, which generates particles with energies uniformly distributed 
between 0 and a specified maximum energy. The particles are initialized at a given altitude with zero 
horizontal velocity and a vertical velocity corresponding to their energy. 
*/

 // UniformSource.hpp


#pragma once
#include "Source.hpp"
#include "constants.h"
#include <random>
#include <sstream>
using namespace std;

class UniformSource : public Source {
public:
  UniformSource() = default;

  void init(const Params& sp,
            vector<float>& z,
            vector<float>& y,
            vector<float>& vz,
            vector<float>& vy,
            vector<float>& E,
            vector<int>&   alive) const override
  {
    // set up RNG
    mt19937                        gen{random_device{}()};
    std::uniform_real_distribution<float> energy_dist(0.0f, sp.Einit);

    // allocate
    z.resize(sp.N);
    y.resize(sp.N);
    vz.resize(sp.N);
    vy.resize(sp.N);
    E.resize(sp.N);
    alive.resize(sp.N);

    for (size_t i = 0; i < sp.N; ++i) {
      float E_sample = energy_dist(gen);

      float v0 = Etov(E_sample); // Convert energy in eV to velocity in m/s

      z[i]     = sp.Zinit;
      y[i]     = 0.0f;
      vz[i]    = -v0;
      vy[i]    = 0.0f;
      E[i]     = E_sample; // in eV
      alive[i] = 1;
    }
  }

  string label(const Params& sp) const override {
    ostringstream os;
    os << "uniform_"
       << static_cast<int>(sp.Einit / 1e3f) << "keV";
    return os.str();
  }
};
