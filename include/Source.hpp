/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the Source class, which is an abstract base class for different types of particle sources 
in the precipitation simulation. The Source class declares a pure virtual method init() that derived classes 
must implement to fill the per-particle arrays with initial conditions. 
It also provides utility methods for energy-velocity conversions and for generating labels and banners
for different source types.
*/

 // Source.hpp

// Source.h
#pragma once
#include <vector>
#include "constants.h"

using namespace constants;

#include "Rnd.h"


struct Params {
  size_t N = 10000;       // Number of primary particles  
  float Zinit;            // Initial position in m
  float Einit;            // energy in eV
  float Emin = 1e3f;      // minimum source sampling energy in eV
  float Emax = 1e6f;      // maximum source sampling energy in eV
  std::string runid = ""; // run identifier, e.g. "TestRun". Default is empty. Used to generate output file names and banners. }; 
};

class Source {
public:
  virtual ~Source() {}
  /// Fill the particle arrays.
  virtual void init(
    const Params& p,
    std::vector<float>& z,      // Vertical position [m]
    std::vector<float>& y,      // Horizontal position [m]
    std::vector<float>& vz,     // Vertical velocity [m/s]
    std::vector<float>& vy,     // Horizontal velocity [m/s]
    std::vector<float>& E,      // Energy [eV]
    std::vector<int>&   alive   // Alive flag (1=alive, 0=dead)
  ) const = 0;

  /// return a filename‐friendly label for this source
  /// e.g. "mono_10keV"
  virtual std::string label(const Params& sp) const = 0;
  /// return a banner string for this source
  virtual std::string banner(const Params& sp) const = 0;


protected:
  float Etov(float E) const {
    // Convert energy in eV to velocity in m/s using relativistic formula
    float gamma = 1.0f + E * constants::eV / (constants::m_e * constants::c * constants::c);
    return constants::c * std::sqrt(1.0f - 1.0f / (gamma * gamma)); // Relativistic
  }

  float vtoE(float v) const {
    // Convert velocity in m/s to energy in eV using relativistic formula
    float gamma = 1.0f / std::sqrt(1.0f - v * v / (constants::c * constants::c));
    return (gamma - 1.0f) * constants::m_e * constants::c * constants::c / constants::eV; // Relativistic
  }

    float EtovNonRel(float E) const {
    // Convert energy in eV to velocity in m/s using non-relativistic formula
    return std::sqrt(2*E*constants::eV/constants::m_e); //  Non-rel
  }

  float vtoENonRel(float v) const {
    // Convert velocity in m/s to energy in eV using non-relativistic formula
    return 0.5f * constants::m_e * v * v / constants::eV;   //Non-rel
  }

};
