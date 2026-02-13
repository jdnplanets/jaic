/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */


 /*
 This file defines the BrokenPowerLawSource class, which implements a source of precipitating particles 
 with a broken power law energy distribution. The source samples particle speeds according to the 
 specified distribution and initializes their state for the precipitation simulation. 
*/
 
 // BrokenPowerLawSource.hpp

#pragma once
#include "Source.hpp"
#include "Rnd.h"
#include "constants.h"
#include <cmath>
#include <sstream>

using namespace std;
using namespace constants;


class BrokenPowerLawSource : public Source {

protected:
  float alpha_, beta_, phi_;

  // helper: sample exactly N speed magnitudes into vvals
  void sampleSpeeds(const Params& sp, vector<float>& vvals) const {

    const float v0 = Etov(phi_); // Convert energy in eV to velocity in m/s

    // Velocity range corresponding to 1–1000 keV
    const float E_min = 1e3f, E_max = 1e6f;
    const float v_min = Etov(E_min);
    const float v_max = Etov(E_max);

    // Target distribution: F(v) ∝ v^3 / [ (v/v0)^α + (v/v0)^β ]
    auto flux_pdf = [&](float v) -> float {
      float ratio = v / v0;
      float denom = pow(ratio, alpha_) + pow(ratio, beta_);
      return (v * v * v) / denom;
    };

    // Estimate max of PDF for rejection sampling
    float fmax = 0.0f;
    for (int i = 0; i < 1000; ++i) {
      float v = v_min + i * (v_max - v_min) / 999;
      fmax = max(fmax, flux_pdf(v));
    }

    // Rejection sampling
    while (vvals.size() < sp.N) {
      float v_try = v_min + rnd() * (v_max - v_min);
      if (rnd() < flux_pdf(v_try) / fmax) {
        vvals.push_back(v_try);
      }
    }
  }

public:
  BrokenPowerLawSource(float alpha, float beta, float phi)
    : alpha_(alpha), beta_(beta), phi_(phi) {}

  void init(const Params& sp,
            vector<float>& z,
            vector<float>& y,
            vector<float>& vz,
            vector<float>& vy,
            vector<float>& E,
            vector<int>&   alive) const override
  {

    vector<float> vvals;
    sampleSpeeds(sp, vvals);

    // Fill particle state
    z.resize(sp.N, sp.Zinit);
    y.resize(sp.N, 0.0f);
    vz.resize(sp.N);
    vy.resize(sp.N, 0.0f);
    E.resize(sp.N);
    alive.resize(sp.N, 1);

    for (size_t i = 0; i < sp.N; ++i){
      vz[i] = -vvals[i];  // downward velocity
      E[i] = vtoE(vvals[i]); // in eV
    }

  }

  string label(const Params& sp) const override {
    ostringstream os;
    os << "BPL_alpha_" << alpha_
       << "_beta_" << beta_
       << "_phi_" << static_cast<int>(phi_) << "V";
    return os.str();
  }

  string banner(const Params& sp) const override {
    ostringstream os;
    os << "Broken Power Law vertical source with alpha = " << alpha_
       << ", beta = " << beta_
       << ", phi = " << phi_ / 1e3f << " kV";
    return os.str();
  }
};
