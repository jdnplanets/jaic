/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines a Gaussian energy distribution source for the precipitation simulation. 
It inherits from the abstract Source class and implements the init() method to fill the per-particle arrays
with initial conditions based on a Gaussian distribution of energies. 
*/


 // GaussianSource.hpp

#pragma once
#include "Source.hpp"
#include "constants.h"
#include <random>
#include <sstream>
using namespace std;

class GaussianSource : public Source {
  float sigma_E_;  // energy‐spread, in eV

public:
  explicit GaussianSource(float sigmaE)
    : sigma_E_(sigmaE)
  {}

  void init(const Params& sp,
            vector<float>& z,
            vector<float>& y,
            vector<float>& vz,
            vector<float>& vy,
            vector<float>& E,
            vector<int>&   alive) const override
  {
    // set up RNG
    mt19937                       gen{random_device{}()};
    normal_distribution<float>    energy_dist(sp.Einit, sigma_E_);

    // allocate
    z.resize(sp.N);
    y.resize(sp.N);
    vz.resize(sp.N);
    vy.resize(sp.N);
    E.resize(sp.N);
    alive.resize(sp.N);

    for (size_t i = 0; i < sp.N; ++i) {
      float E_sample = energy_dist(gen);
      // guard against negative draws
      if (E_sample < 0.0f) E_sample = 0.0f;
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
    os << "gaussE_";
    if (sp.Einit >= 1e3)
      os << static_cast<int>(sp.Einit / 1e3) << "keV";
    else
      os << static_cast<int>(sp.Einit) << "eV";

    os << "_sigma_";
    if (sigma_E_ >= 1e3)
      os << static_cast<int>(sigma_E_ / 1e3) << "keV";
    else
      os << static_cast<int>(sigma_E_) << "eV";

    return os.str();
  }

std::string banner(const Params& sp) const override {
    ostringstream os;
    os << "Gaussian source with <E> = " << sp.Einit/1e3 << " keV"
       << " and sigma = " << sigma_E_ << " eV" << endl;

    return os.str();
  }
};
