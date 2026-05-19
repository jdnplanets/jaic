/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
World1D allocates the grids/arrays (based on geometry), but does not decide how to populate them.
Derived classes override populate() to fill nH2/T/hydrocarbons/etc (from files, analytic, etc.).
*/

// World1D.hpp

#pragma once

#include "ArrayOps.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>

// Generic override bag
#include "WorldOverrides.hpp"

struct IIonChemistryModel;  // forward decl

struct World1D {

    // Output directory for generated data files
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    static inline std::string datadir = jaicroot + "data/world1d/";

    // -----------------------------
    // Geometry / grid
    // -----------------------------
    float P0 = 1e5f; // Reference altitude pressure in Pa (1 bar)
    float P1 = 1e-7f; // Top of atmosphere pressure in Pa (1e-12 bar is ~3000 km for Jupiter)
    int nP = 500; // Number of pressure bins
    std::vector<float> P; // Pressure grid in Pa, size nP


    float Z1;  // top of domain (m)
    float Z0;     // bottom of domain (m)
    int   nz = nP;      // number of altitude bins

    std::vector<float> Z;  // altitude grid (m), size nz
    std::vector<float> dz;     // bin size (m)   computed
    std::vector<float> dzcm;     // bin size (cm)  computed
    std::vector<float> Z_edges; // altitude edges (m), size nz+1 computed

    // -----------------------------
    // Geometry of input files (precomputed)
    // -----------------------------
    static inline const int nPin = 500; // number of pressure bins in the input atmosphere file (atmosphere.dat)
    static inline const float P0_in = 1e5f;   
    static inline const float P1_in = 1e-12f;
    std::vector<float> P_in; // pressure grid of input file (bar)

    // -----------------------------
    // State vectors (size nP)
    // -----------------------------
    // Neutral densities (cm^-3)
    std::vector<float> nH2; // read in as floats for GPU compatability
    std::vector<float> nHe;
    std::vector<float> nH;
    std::vector<float> nCH4;
    std::vector<float> nC2H2;
    std::vector<float> nC2H4;
    std::vector<float> nC2H6;

    // Ion & electron densities (cm^-3)
    std::vector<float> ne;
    std::vector<float> nH3p;
    std::vector<float> nCH5p;
    std::vector<float> nC3Hnp;

    // Temperature (K)
    std::vector<float> T;


    // Magnetic field (Tesla)
    double B = 0.001;
    double Bdipang = 90.0; // dip angle of the magnetic field in degrees (90 = vertical, 0 = horizontal)

    // Ion chemistry model (chosen per-world in WorldFactory).
    // If null, IonDensity can fall back to a default (e.g. all-H3+).
    std::shared_ptr<const IIonChemistryModel> ionChemistryModel;
    // Destructor
    virtual ~World1D() = default;

    // Base constructor: set geometry, build Z, allocate arrays.
    // Derived classes "override" Z0/Z1/nz by calling this ctor with different values.
    explicit World1D(float P0_bar = 1.0f, float P1_bar = 1e-12f, int nP_ = 500)
        : P0(P0_bar), P1(P1_bar), nP(nP_)
    {
        if (nP <= 0) {
            throw std::invalid_argument("World1D: nP must be > 0");
        }

        // Allocate arrays (populate() will fill values)
        P.assign(nP, 0.0f);
        Z.assign(nP, 0.0f);
        dz.assign(nP, 0.0f);
        dzcm.assign(nP, 0.0f);
        Z_edges.assign(nP + 1, 0.0f);

        nH2.assign(nP, 0.0f);
        nHe.assign(nP, 0.0f);
        nH.assign(nP, 0.0f);
        nCH4.assign(nP, 0.0);
        nC2H2.assign(nP, 0.0);
        nC2H4.assign(nP, 0.0);
        nC2H6.assign(nP, 0.0);

        ne.assign(nP, 0.0f);
        nH3p.assign(nP, 0.0f);
        nCH5p.assign(nP, 0.0f);
        nC3Hnp.assign(nP, 0.0f);

        T.assign(nP, 0.0);
    }

    // Apply optional overrides init(). Base implementation applies only generic knobs.
    // Derived classes override this to apply e.g. Tiso_K, g, n0, etc.
    virtual void apply_overrides(const WorldOverrides& o) {
        if (o.B_T) {
            B = *o.B_T;
        }
        if (o.Bdipang_deg) {
            Bdipang = *o.Bdipang_deg;
        }
    }

    // Two-stage init: call once after construction (e.g. in factory or main).
    void init() {
        populate();
        validate();
        std::cout << "World initialized with " << nP << " pressure levels from " << P0 << " to " << P1 << " Pa." << std::endl;
        std::cout << "Altitude range: " << Z0/1e3 << " km to " << Z1/1e3 << " km." << std::endl;

        // Copy all lower edges
        Z_edges = Z;
         // Infer the final upper edge from the last spacing
        const float dz_top = Z.back() - Z[nz - 2];
        Z_edges.push_back(Z.back() + dz_top);
        
    }


    void write_species_densities(const std::string& output_file) const {
        const std::string fullpath = (std::filesystem::path(datadir) / output_file).string();
        std::ofstream file(fullpath);
        if (!file) {
            throw std::runtime_error("Jupiter: failed to open output file: " + fullpath);
        }

        // Write header
        file << std::left
            
             << std::setw(12) << "z_m"
             << std::setw(12) << "nH2"
             << std::setw(12) << "nCH4"
             << std::setw(12) << "nC2H2"
             << std::setw(12) << "nC2H4"
             << std::setw(12) << "nC2H6"
             << std::setw(12) << "T_K" << "\n";


        // Write data rows
        for (int i = 0; i < nP; ++i) {
            file << std::left
                <<  std::setw(12) << P[i]
                 << std::setw(12) << Z[i]
                 << std::setw(12) << nH2[i]
                 << std::setw(12) << nHe[i]
                 << std::setw(12) << nH[i]
                 << std::setw(12) << nCH4[i]
                 << std::setw(12) << nC2H2[i]
                 << std::setw(12) << nC2H4[i]
                 << std::setw(12) << nC2H6[i]
                 << std::setw(12) << T[i] << "\n";
        }
    }

protected:
    // Derived world must populate arrays (file I/O or analytic).
    virtual void populate() = 0;

    // Optional sanity checks
    virtual void validate() const {
        if (static_cast<int>(Z.size()) != nP)   throw std::runtime_error("World1D: Z wrong size");
        if (static_cast<int>(nH2.size()) != nP) throw std::runtime_error("World1D: nH2 wrong size");
        if (static_cast<int>(nHe.size()) != nP) throw std::runtime_error("World1D: nHe wrong size");
        if (static_cast<int>(nH.size()) != nP)  throw std::runtime_error("World1D: nH wrong size");
        if (static_cast<int>(nCH4.size()) != nP) throw std::runtime_error("World1D: nCH4 wrong size");
        if (static_cast<int>(nC2H2.size()) != nP) throw std::runtime_error("World1D: nC2H2 wrong size");
        if (static_cast<int>(nC2H4.size()) != nP) throw std::runtime_error("World1D: nC2H4 wrong size");
        if (static_cast<int>(nC2H6.size()) != nP) throw std::runtime_error("World1D: nC2H6 wrong size");
        if (static_cast<int>(ne.size()) != nP)   throw std::runtime_error("World1D: ne wrong size");
        if (static_cast<int>(nH3p.size()) != nP)  throw std::runtime_error("World1D: nH3p wrong size");
        if (static_cast<int>(nCH5p.size()) != nP)  throw std::runtime_error("World1D: nCH5p wrong size");
        if (static_cast<int>(nC3Hnp.size()) != nP)  throw std::runtime_error("World1D: nC3Hnp wrong size");
        if (static_cast<int>(T.size()) != nP)   throw std::runtime_error("World1D: T wrong size");
    }

    std::string path_join(const std::string& dir, const std::string& file) const {
        if (dir.empty()) return file;
        if (dir.back() == '/' || dir.back() == '\\') return dir + file;
        return dir + "/" + file;
    }
};


