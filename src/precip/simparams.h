/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the SimParams struct, which contains all the parameters for the simulation, 
including the energy grid definition and convenience functions for working with the energy grid. 
It also includes parameters for the precipitation simulation, such as the number of particles, 
initial energy, and initial position. The energy grid can be defined as either logarithmic or linear, 
and the SimParams struct provides functions to convert between energies and their corresponding indices 
in the energy grid.

*/

 // simparams.h

#pragma once

#include "Egrid.hpp"
#include <algorithm> // for std::min/std::max (host)

struct SimParams {

    float P0;           // bottom of domain, Pa
    float P1;           // top of domain, Pa
    int nbinsP;         // Number of pressure bins. These are set by World1d
    float Z1;           // Top of the simulation domain in m (defines distance units)
    float Z0;           // Bottom of the simulation domain in m
    int nbinsZ;         // Number of altitude bins, derived from Z0, Z1, and dz. These are set by World1d
    float dz;           // Size of each bin in m, derived from Z0, Z1, and nbinsZ. 
    float dzcm;         // in cm

    // Used energy grid definition
    float E0      = 1.0f;                         // bottom of energy range in eV (edge 0)
    float Emax;                                   // top of energy range in eV (edge nbinsE); derived from E0/decades
    static constexpr float decades = 6.0f;        // number of decades in energy range (defines Emax = E0*10^decades)
    static constexpr int nE        = 500;         // Number of energy bins in the energy grid
    EgridType egridType = EgridType::Logarithmic; // type of energy grid (logarithmic or linear)

    //Energy grid of input files (precomputed)
    static constexpr float E0_in      = 1.0f;        // bottom of energy range in eV (edge 0)
    static constexpr float Emax_in    = 1e6f;        // top of energy range in eV (edge nbinsE); derived from E0/decades
    static constexpr float decades_in = 6.0f;        // number of decades in energy range (defines Emax = E0*10^decades)
    static constexpr int nE_in        = 100;         // Number of energy bins in the energy grid
    EgridType egridType_in = EgridType::Logarithmic; // type of energy grid (logarithmic or linear) 


    // Angular / energy discretization
    static constexpr int nTheta = 1000;  // Number of bins in scattering angle PDF (theta)
    static constexpr int nbinsE = nE;    // Alias


    // Precipitation parameters
    static constexpr int nSpecies = 1;  // Number of scattering species
    int          N;                     // Number of primary particles, populated from the source
    float        Einit;                 // Initial energy in eV (only used for monoenergetic sources but required anyway)
    float        Zinit;                 // Initial position in m
    float        maxdt        = 1e-5f;  // Maximum time step in seconds used for the simulation
    const float  Emin         = 7.0f;   // Minimum energy of alive electrons in eV


    const float  cosThetaMean = 0.64f;   // Mean of the cosine of the scattering angle for the two-stream approximation


    // Collision typees included
    enum Collision {
        elastic,
        ionisation,
        excitation_triplet_a,
        excitation_triplet_b,
        excitation_triplet_c,
        excitation_triplet_e,
        excitation_singlet_B,
        excitation_singlet_C,
        vibrational,
        rotational,
        excitation_singlet_EF,
        collisionCount // automatically equals the number of collisions
    };
    static constexpr int nCollisions = collisionCount;

    // Energy grid convenience functions
    
    EGRID_HD inline void update_derived() {
        if (E0 <= 0.0f) E0 = 1.0f;
        Emax = E0 * powf(10.0f, decades);
        if (Emax <= E0) Emax = E0 * 10.0f; // safety
    }

    // Convenience: the extra parameter needed by typed statics
    // (decades for log grids; Emax for linear grids)
    EGRID_HD inline float EgridParam() const {
        return (egridType == EgridType::Logarithmic) ? decades : Emax;
    }

    // Small clamps to keep callers safe
    EGRID_HD inline int clamp_edge(int edge) const {
        if (edge < 0) return 0;
        if (edge > nbinsE) return nbinsE;
        return edge;
    }
    EGRID_HD inline int clamp_bin(int e) const {
        if (e < 0) return 0;
        if (e >= nbinsE) return nbinsE - 1;
        return e;
    }

    // ---------- Indexing ----------
    EGRID_HD inline int EToe(float E) const {
        return Egrid::EToIndex_static_typed(E, nbinsE, E0, EgridParam(), egridType);
    }

    EGRID_HD inline int Etoe(float E) const {
        return Egrid::EToIndex_static_typed(E, nbinsE, E0, EgridParam(), egridType);
    }

    // ---------- Edge energies ----------
    EGRID_HD inline float El(int edge) const {   // edge in [0, nbinsE]
        edge = clamp_edge(edge);
        return Egrid::eToE_static_typed(edge, nbinsE, E0, EgridParam(), egridType);
    }

    // ---------- Bin centres ----------
    EGRID_HD inline float Ec(int e) const {        // "default centre": Log->geom, Linear->arith
        e = clamp_bin(e);
        return Egrid::eToEc_static_typed(e, nbinsE, E0, EgridParam(), egridType);
    }

    EGRID_HD inline float Ecg(int e) const {        // geometric mean of edges
        e = clamp_bin(e);
        return Egrid::eToEcg_static_typed(e, nbinsE, E0, EgridParam(), egridType);
    }

    EGRID_HD inline float Eca(int e) const {        // arithmetic mean of edges
        e = clamp_bin(e);
        return Egrid::eToEca_static_typed(e, nbinsE, E0, EgridParam(), egridType);
    }

    // ---------- Bin width ----------
    EGRID_HD inline float dE(int e) const {
        e = clamp_bin(e);
        return Egrid::eTodE_static_typed(e, nbinsE, E0, EgridParam(), egridType);
    }

    // Altitude grid convenience functions (reuse EGRID_HD for host/device compatibility))

    // Must be set to point to host memory when called on host,
    // and device memory when called on device.
    const float* Z_edges = nullptr;  // size nbinsZ + 1

    EGRID_HD inline int Ztoz(float Z) const {

        if (Z_edges == nullptr || nbinsZ <= 0) return 0;

        if (Z <= Z_edges[0]) return 0;
        if (Z >= Z_edges[nbinsZ]) return nbinsZ - 1;

        // Binary search to find the right bin
        int lo = 0;
        int hi = nbinsZ;

        

        while (hi - lo > 1) {
            int mid = lo + (hi - lo) / 2;

            if (Z_edges[mid] <= Z) {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        return lo;
    }


};
