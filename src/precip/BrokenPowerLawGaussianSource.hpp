#pragma once
#include "Source.hpp"
#include "Rnd.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

using namespace std;
using namespace constants;


class BrokenPowerLawGaussianSource : public Source {

protected:
  float alpha_, beta_, phi_;
  float sigma_;          // Gaussian sigma in eV
  float peakFactor_;     // Total at E=phi_ is this factor above pure BPL (default 5)


  float sigmaV() const {
    const float v0_ = Etov(phi_);
    if (sigma_ <= 0.0f || v0_ <= 0.0f) return 0.0f;
    // sigma_E [eV] -> sigma_v [m/s]
    return sigma_ * eV / (m_e * v0_);
  }
  
  // Physical distribution function f(v) up to an overall constant f0
  float physicalFv(float v) const {
    const float v0_ = Etov(phi_);
    const float ratio = v / v0_;
    const float bpl = 1.0f / (pow(ratio, alpha_) + pow(ratio, beta_));

    if (sigma_ <= 0.0f || peakFactor_ <= 1.0f) {
      return bpl;
    }

    const float sigv = sigmaV();
    const float x = (v - v0_) / sigv;

    // At v=v0, f_BPL(v0)=1/2, so Gaussian amplitude must be
    // (peakFactor-1)/2 to make total = peakFactor * f_BPL(v0)
    const float gaussAmp = 0.5f * (peakFactor_ - 1.0f);
    const float gauss = gaussAmp * exp(-0.5f * x * x);

    return bpl + gauss;
  }

    // Sampling PDF in speed magnitude: p(v) ∝ v^3 f(v)
  float samplingPdf(float v) const {
    return v * v * v * physicalFv(v);
  }

  // helper: sample exactly N speed magnitudes into vvals
  void sampleSpeeds(const Params& sp, vector<float>& vvals) const {

    // Velocity range from source parameters.
    const float E_min = sp.Emin;
    const float E_max = std::max(sp.Emax, E_min + 1.0f);
    const float v_min = Etov(E_min);
    const float v_max = Etov(E_max);

    // Estimate max of composite PDF for rejection sampling
    float fmax = 0.0f;
    for (int i = 0; i < 4000; ++i) {
      const float v = v_min + i * (v_max - v_min) / 3999.0f;
      fmax = max(fmax, samplingPdf(v));
    }

    // Safety fallback
    if (fmax <= 0.0f) fmax = 1.0f;

    // Rejection sampling
    vvals.reserve(sp.N);
    while (vvals.size() < static_cast<size_t>(sp.N)) {
      const float v_try = v_min + rnd() * (v_max - v_min);
      if (rnd() < samplingPdf(v_try) / fmax) {
        vvals.push_back(v_try);
      }
    }
  }

public:
  BrokenPowerLawGaussianSource(float alpha,
                               float beta,
                               float phi,
                               float sigma,
                               float peakFactor = 5.0f)
    : alpha_(alpha),
      beta_(beta),
      phi_(phi),
      sigma_(sigma),
      peakFactor_(peakFactor) {}

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

    z.resize(sp.N, sp.Zinit);
    y.resize(sp.N, 0.0f);
    vz.resize(sp.N);
    vy.resize(sp.N, 0.0f);
    E.resize(sp.N);
    alive.resize(sp.N, 1);

    for (size_t i = 0; i < static_cast<size_t>(sp.N); ++i) {
      vz[i] = -vvals[i];   // downward velocity
      E[i] = vtoE(vvals[i]); // in eV
    }
  }

  string label(const Params& sp) const override {
    ostringstream os;
    os << "BPLG_alpha_" << alpha_
       << "_beta_" << beta_
       << "_phi_" << static_cast<int>(phi_) << "V"
       << "_sigma_" << static_cast<int>(sigma_) << "V"
       << "_peak_" << peakFactor_;
    return os.str();
  }

  string banner(const Params& sp) const override {
    ostringstream os;
    os << "Broken Power Law + Gaussian vertical source with alpha = " << alpha_
       << ", beta = " << beta_
       << ", phi = " << phi_ / 1e3f << " kV"
       << ", sigma = " << sigma_ / 1e3f << " keV"
       << ", peak factor = " << peakFactor_;
    return os.str();
  }
};