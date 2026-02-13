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

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>
#include <filesystem>

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
    float Z1   = 3000e3f;  // top of domain (m)
    float Z0   = 0.0f;     // bottom of domain (m)
    int   nz   = 300;      // number of altitude bins

    float dz   = 0.0f;     // bin size (m)   computed
    float dzcm = 0.0f;     // bin size (cm)  computed

    std::vector<float> Z;  // altitude grid (m), size nz

    // -----------------------------
    // State vectors (size nz)
    // -----------------------------
    // Neutral densities (cm^-3)
    std::vector<float>  nH2;
    std::vector<double> nCH4;
    std::vector<double> nC2H2;
    std::vector<double> nC2H4;
    std::vector<double> nC2H6;

    // Ion & electron densities (cm^-3)
    std::vector<float> ne;
    std::vector<float> nH3p;
    std::vector<float> nCH5p;
    std::vector<float> nC3Hnp;

    // Temperature (K)
    std::vector<double> T;

    // Magnetic field (Tesla)
    double B = 0.0015;


    // Ion chemistry model (chosen per-world in WorldFactory).
    // If null, IonDensity can fall back to a default (e.g. all-H3+).
    std::shared_ptr<const IIonChemistryModel> ionChemistryModel;
    // Destructor
    virtual ~World1D() = default;

    // Base constructor: set geometry, build Z, allocate arrays.
    // Derived classes "override" Z0/Z1/nz by calling this ctor with different values.
    explicit World1D(float Z0_m = 0.0f, float Z1_m = 3000e3f, int nz_ = 300)
        : Z1(Z1_m), Z0(Z0_m), nz(nz_)
    {
        if (nz <= 0) {
            throw std::invalid_argument("World1D: nz must be > 0");
        }

        dz   = (Z1 - Z0) / static_cast<float>(nz);
        dzcm = dz * 1e2f;

        // Build Z grid
        Z.resize(nz);
        for (int i = 0; i < nz; ++i) {
            Z[i] = Z0 + static_cast<float>(i) * dz;
        }

        // Allocate arrays (populate() will fill values)
        nH2.assign(nz, 0.0f);
        nCH4.assign(nz, 0.0);
        nC2H2.assign(nz, 0.0);
        nC2H4.assign(nz, 0.0);
        nC2H6.assign(nz, 0.0);

        ne.assign(nz, 0.0f);
        nH3p.assign(nz, 0.0f);
        nCH5p.assign(nz, 0.0f);
        nC3Hnp.assign(nz, 0.0f);

        T.assign(nz, 0.0);
    }

    // Apply optional overrides init(). Base implementation applies only generic knobs.
    // Derived classes override this to apply e.g. Tiso_K, g, n0, etc.
    virtual void apply_overrides(const WorldOverrides& o) {
        if (o.B_T) {
            B = *o.B_T;
        }
    }

    // Two-stage init: call once after construction (e.g. in factory or main).
    void init() {
        populate();
        validate();
    }




    void write_hydrocarbon_densities(const std::string& output_file) const {
        const std::string fullpath = (std::filesystem::path(datadir) / output_file).string();
        std::ofstream file(fullpath);
        if (!file) {
            throw std::runtime_error("Jupiter: failed to open output file: " + fullpath);
        }

        // Write header
        file << "z_m,nCH4,nC2H2,nC2H4,nC2H6\n";

        // Write data rows
        for (int i = 0; i < nz; ++i) {
            file << Z[i] << "," 
                 << nCH4[i] << "," 
                 << nC2H2[i] << "," 
                 << nC2H4[i] << "," 
                 << nC2H6[i] << "\n";
        }
    }

protected:
    // Derived world must populate arrays (file I/O or analytic).
    virtual void populate() = 0;

    // Optional sanity checks
    virtual void validate() const {
        if (static_cast<int>(Z.size()) != nz)   throw std::runtime_error("World1D: Z wrong size");
        if (static_cast<int>(nH2.size()) != nz) throw std::runtime_error("World1D: nH2 wrong size");
        if (static_cast<int>(T.size()) != nz)   throw std::runtime_error("World1D: T wrong size");
    }

    std::string path_join(const std::string& dir, const std::string& file) const {
        if (dir.empty()) return file;
        if (dir.back() == '/' || dir.back() == '\\') return dir + file;
        return dir + "/" + file;
    }
};


