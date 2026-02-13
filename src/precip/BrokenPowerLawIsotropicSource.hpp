/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the BrokenPowerLawIsotropicSource class, which is a specific implementation 
of the BrokenPowerLawSource. It generates an isotropic distribution of particle velocities 
following a broken power law in speed, with parameters alpha, beta, and phi.
 */

 // BrokenPowerLawIsotropicSource.hpp


#pragma once
#include "BrokenPowerLawSource.hpp"
#include "constants.h"
#include "Rnd.h"
#include <cmath>
#include <sstream>

using namespace std;
using namespace constants;


class BrokenPowerLawIsotropicSource : public BrokenPowerLawSource {
  float alpha_, beta_, phi_;

public:
  BrokenPowerLawIsotropicSource(float alpha, float beta, float phi)
    : BrokenPowerLawSource(alpha, beta, phi),   // Initialize the base class
      alpha_(alpha), beta_(beta), phi_(phi) {}

  void init(const Params& sp,
            vector<float>& z,
            vector<float>& y,
            vector<float>& vz,
            vector<float>& vy,
            vector<float>& E,
            vector<int>&   alive) const override
  {
    // get N speed magnitudes
    vector<float> vvals;
    sampleSpeeds(sp, vvals);

    // Fill particle state
    z.resize(sp.N, sp.Zinit);  // in m
    y.resize(sp.N, 0.0f);
    vz.resize(sp.N);           // in m/s
    vy.resize(sp.N, 0.0f);
    E.resize(sp.N, 0.0f);      // in eV
    alive.resize(sp.N, 1);

    for (size_t i = 0; i < sp.N; ++i) {
        float theta = pi * rnd() - pi/2; // uniform distribution
        vz[i] = -vvals[i] * cos(theta);
        vy[i] = vvals[i] * sin(theta);
        E[i] = vtoE(vvals[i]); // set the energy in eV
    }

  }
  

  string label(const Params& sp) const override {
    ostringstream os;
    os << "BPLiso_alpha_" << alpha_
       << "_beta_" << beta_
       << "_phi_" << static_cast<int>(phi_) << "V";
    return os.str();
  }

  string banner(const Params& sp) const override {
    ostringstream os;
    os << "Broken Power Law isotropic source with alpha = " << alpha_
       << ", beta = " << beta_
       << ", phi = " << phi_ / 1e3f << " kV";
    return os.str();
  }

};