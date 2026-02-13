/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the DataSource class, which implements the Source interface to provide a particle source 
based on input data from a file. The file should contain two columns: energy (in keV) and flux (in cm^-2 s^-1 sr^-1 keV^-1). 
The DataSource class reads this data, constructs a log-spaced energy grid, and computes weights for sampling 
particle energies according to the provided flux distribution. It also accounts for an isotropic pitch angle 
distribution within a specified half-width cone.
*/

 // DataSource.hpp


#pragma once
#include "Source.hpp"
#include "constants.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>

using namespace std;
using namespace constants;



class DataSource : public Source {
  vector<float>   Egrid_;     // keV bin edges, size = Ngrid+1
  vector<double>  weights_;   // weights per bin, size = Ngrid
  float           thetaMax_;  // pitch‐angle half‐width (rad)
  string          filename_;  // store filename

private:
  double nFlux = 0.0; // Number flux in cm-2 s-1
  double EFlux = 0.0; // Energy flux in mW m-2

public:
  /// filename: two cols [E_keV, flux_cm2_s_sr_keV]
  /// thetaHalfDeg: half‐width in degrees (default 90°)
  /// Ngrid: number of log‐energy bins
  DataSource(const string& filename,
                     float thetaHalfDeg = 90.0f,
                     size_t Ngrid = 200)
    : thetaMax_(thetaHalfDeg * pi/180.0f), filename_(filename)
  {
    // Read raw data
    static  const string home = getenv("HOME");
    static  const string jaicroot = home + "/cpp/jaic/";
    static  const string outdir = jaicroot + "input/";
    ifstream in(outdir+filename);
    if (!in) throw runtime_error("Cannot open " + filename);
    cout << "Reading input file: " << outdir + filename << endl;
    vector<float> Eraw, Fraw;
    for (float E, F; in >> E >> F; ) {
      if (E > 0 && F > 0) {
      Eraw.push_back(E);
      Fraw.push_back(F);
      }
    }
    if (Eraw.size() < 2) throw runtime_error("Not enough data");

    // Sort by E
    vector<size_t> idx(Eraw.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
      [&](size_t i, size_t j){ return Eraw[i]<Eraw[j]; });
    vector<float> Esorted, Fsorted;
    for (auto i: idx) { Esorted.push_back(Eraw[i]); Fsorted.push_back(Fraw[i]); }

    // Build log‐grid
    float Emin = Esorted.front(), Emax = Esorted.back();
    Egrid_.resize(Ngrid+1);
    float logEmin = log(Emin), logEmax = log(Emax);
    for (size_t i=0; i<=Ngrid; ++i)
      Egrid_[i] = exp(logEmin + (logEmax-logEmin)*i/Ngrid);

    // Interpolate F onto grid mid-points & build weights
    weights_.resize(Ngrid);
    // solid angle of a cone: 2π(1–cosθ)
    float solidAngle = pi*(1.0f - cos(thetaMax_)*cos(thetaMax_));

    for (size_t i=0; i<Ngrid; ++i) {
      float E0 = Egrid_[i], E1 = Egrid_[i+1];
      float Emid = sqrt(E0*E1);
      // locate in sorted data
      auto it = lower_bound(Esorted.begin(), Esorted.end(), Emid);
      double Fmid;
      if (it==Esorted.begin()) {
        Fmid = Fsorted.front();
      } else if (it==Esorted.end()) {
        Fmid = Fsorted.back();
      } else {
        size_t j = (it - Esorted.begin());
        // log-log interpolation
        float E_lo = Esorted[j-1], E_hi = Esorted[j];
        float F_lo = Fsorted[j-1], F_hi = Fsorted[j];
        float t = (log(Emid)-log(E_lo)) / (log(E_hi)-log(E_lo));
        Fmid = exp(log(F_lo) + t*(log(F_hi)-log(F_lo)));
      }
      float dE = E1 - E0;
      // weight ∝ flux * dE * solid angle
      weights_[i] = Fmid * dE * solidAngle;
    }

    // ormalize weights to sum = 1
    vector<double> unnormed_weights_ = weights_;
    double sumW = accumulate(weights_.begin(), weights_.end(), 0.0);
    for (auto& w : weights_) w /= sumW;
    
    nFlux = sumW;
    cout << "Total number flux (multiply Precip results by this): " << nFlux << " cm-2 s-1" << endl;

    for (size_t i = 0; i < Ngrid; ++i) {
      float E0 = Egrid_[i], E1 = Egrid_[i+1];
      float Emid = sqrt(E0 * E1);
      EFlux += unnormed_weights_[i] * Emid;
    }
    EFlux *= 1e3 * qe * 1e4 * 1e3; // in mW m-2
    cout << "DataSource: Total energy flux = " << EFlux << " mW m-2" << endl;

  }

  void init(const Params& sp,
            vector<float>& z,
            vector<float>& y,
            vector<float>& vz,
            vector<float>& vy,
            vector<float>& E,
            vector<int>&   alive) const override
  {
    // RNG
    static thread_local mt19937 gen{random_device{}()};
    discrete_distribution<size_t> dist(weights_.begin(), weights_.end());

    // Prepare output
    z.resize(sp.N, sp.Zinit);
    y.resize(sp.N, 0.0f);
    vz.resize(sp.N);
    vy.resize(sp.N, 0.0f);
    E.resize(sp.N);
    alive.resize(sp.N, 1);

    for (size_t i = 0; i < sp.N; ++i) {
      // --- sample energy bin ---
      size_t bin = dist(gen);
      float E0 = Egrid_[bin], E1 = Egrid_[bin+1];
      // sample log‐uniform within [E0,E1]
      float t = rnd(); // uniform [0,1)
      float E_keV = E0 * pow(E1/E0, t);

      // convert to speed
      float E_eV  = E_keV * 1e3f;
      E[i]       = E_eV; // in eV
      float v    = Etov(E_eV); // in m/s

      // --- sample pitch angle within cone ---
      float theta = thetaMax_ * rnd() - thetaMax_/2; // uniform distribution
      vz[i] = -v * cos(theta);
      vy[i] =  v * sin(theta);

    }

  }

  std::string label(const Params& sp) const override {
    ostringstream os;
    os << "fileSpecθ" << (thetaMax_*180.0f/pi) << "deg";
    return os.str();
  }

  std::string banner(const Params& sp) const override {
    ostringstream os;
    os << "Isotropic Data source from file = " << filename_
       << " with θ_half = " << (thetaMax_*180.0f/pi) << "°" << endl
       << "  ----> Number flux = " << scientific << nFlux << " cm-2 s-1 (the value by which to multiply Precip results)" << endl
       << "  ----> Energy flux = " << defaultfloat << EFlux << " mW m-2" << endl
       << "Use src->getNumberFlux() to get the number flux in cm-2 s-1" << endl;

    return os.str();
  }

  double getNumberFlux() const {
    if (nFlux <= 0.0) {
      throw runtime_error("Number flux is not set or is zero.");
    }
    return nFlux; // in cm-2 s-1
  }


};
