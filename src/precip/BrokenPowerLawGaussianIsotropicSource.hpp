// FluxWeightedSource.h
#pragma once
#include "BrokenPowerLawGaussianSource.hpp"
#include "constants.h"
#include "Rnd.h"
#include <cmath>
#include <sstream>

using namespace std;
using namespace constants;


class BrokenPowerLawGaussianIsotropicSource : public BrokenPowerLawGaussianSource {
  float alpha_, beta_, phi_, sigma_, peakFactor_;

public:
  BrokenPowerLawGaussianIsotropicSource(float alpha, float beta, float phi, float sigma, float peakFactor = 5.0f)
    : BrokenPowerLawGaussianSource(alpha, beta, phi, sigma, peakFactor),   // Initialize the base class
      alpha_(alpha), beta_(beta), phi_(phi), sigma_(sigma), peakFactor_(peakFactor) {}

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
    z.resize(sp.N, sp.Zinit);
    y.resize(sp.N, 0.0f);
    vz.resize(sp.N);
    vy.resize(sp.N, 0.0f);
    E.resize(sp.N);
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
    os << "BPLGiso_alpha_" << alpha_
       << "_beta_" << beta_
       << "_phi_" << static_cast<int>(phi_) << "V"
        << "_sigma_" << static_cast<int>(sigma_) << "V"
        << "_peak_" << peakFactor_;
    return os.str();
  }

  string banner(const Params& sp) const override {
    ostringstream os;
    os << "Broken Power Law + Gaussian isotropic source with alpha = " << alpha_
       << ", beta = " << beta_
       << ", phi = " << phi_ / 1e3f << " kV"
       << ", sigma = " << sigma_ / 1e3f << " keV"
       << ", peak factor = " << peakFactor_;
    return os.str();
  }

};