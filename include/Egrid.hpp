/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // Egrid.hpp

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

#ifndef EGRID_HPP
#define EGRID_HPP

// Host/device annotation for NVCC
#if defined(__CUDACC__)
#define EGRID_HD __host__ __device__
#else
#define EGRID_HD
#endif

enum class EgridType { Logarithmic, Linear };

// Forward declarations for concrete helper structs
struct LogEgrid;
struct LinEgrid;

// Base Egrid: lightweight value type that dispatches to concrete helpers.
//   - nbins is the number of energy bins
//   - There are nbins+1 edges (grid points).
//   - eToE(e) returns the EDGE energy for edge index e in [0, nbins].
//   - eToEca/eToEcg/eTodE take bin index e in [0, nbins-1].
//
class Egrid {
public:
    int nbins;      // number of energy bins
    double decades; // for log grids: log10(Emax/Emin)
    double Emin;    // minimum energy (eV) = edge 0
    double Emax;    // maximum energy (eV) = edge nbins
    EgridType type; // grid flavour


    std::string typestring() const {
        switch (type) {
            case EgridType::Logarithmic: return "logarithmic";
            case EgridType::Linear: return "linear";
            default: return "unknown";
        }
    }

    Egrid(int nbins_ = 100, double decades_ = 6.0, double Emin_ = 1.0, EgridType t = EgridType::Logarithmic)
        : nbins(std::max(2, nbins_)), decades(decades_), Emin(Emin_), Emax(0.0), type(t) {
        if (Emin <= 0.0) Emin = 1.0;
        Emax = Emin * std::pow(10.0, decades);
    }

    // Factory from explicit Emin/Emax
    static inline Egrid fromEminEmax(int nbins_, double Emin_, double Emax_, EgridType t = EgridType::Logarithmic) {
        if (nbins_ < 2) nbins_ = 2;
        if (Emin_ <= 0.0) Emin_ = 1.0;
        if (Emax_ <= Emin_) Emax_ = Emin_ * 1e6;
        double decades_ = std::log10(Emax_ / Emin_);
        Egrid g(nbins_, decades_, Emin_, t);
        g.Emax = Emax_;
        return g;
    }

    // Edge energy at index e in [0, nbins] (inclusive)
    inline double eToE(int e) const;
    // geometric mean (center) of BIN e in [0, nbins-1]
    inline double eToEcg(int e) const;
    // arithmetic mean (center) of BIN e in [0, nbins-1]
    inline double eToEca(int e) const;
    // backward-compatible alias
    inline double eToEc(int e) const;
    // bin width of BIN e in [0, nbins-1]
    inline double eTodE(int e) const;
    // map energy to BIN index in [0, nbins-1]
    inline int EToIndex(double E) const;
   // Return all bin edges: LH edge of every bin plus the final RH edge.
    // Size is nbins + 1; edges[i] == eToE(i), for i=0..nbins.
    inline std::vector<double> getEdges() const; 

    // Backwards-compatible static alias: default is logarithmic
    EGRID_HD static inline int EToIndex_static(double E, int nbins_, double Emin_, double decades_);
    EGRID_HD static inline int EToIndex_static(float E, int nbins_, float Emin_, float decades_);


    // Explicit linear static helpers
    EGRID_HD static inline int EToIndex_static_linear(double E, int nbins_, double Emin_, double Emax_);
    EGRID_HD static inline int EToIndex_static_linear(float  E, int nbins_, float  Emin_, float  Emax_);

    // Typed static dispatch: choose mapping by grid type
    EGRID_HD static inline int EToIndex_static_typed(double E, int nbins_, double Emin_, double decades_or_Emax_,
                                                     EgridType t);
    EGRID_HD static inline int EToIndex_static_typed(float  E, int nbins_, float  Emin_, float  decades_or_Emax_,
                                                     EgridType t);

        // --- Static helpers (logarithmic) ---
    EGRID_HD static inline double eToE_static_log  (int e, int nbins_, double Emin_, double decades_);
    EGRID_HD static inline double eToEc_static_log (int e, int nbins_, double Emin_, double decades_);
    EGRID_HD static inline double eToEcg_static_log(int e, int nbins_, double Emin_, double decades_);
    EGRID_HD static inline double eToEca_static_log(int e, int nbins_, double Emin_, double decades_);
    EGRID_HD static inline double eTodE_static_log (int e, int nbins_, double Emin_, double decades_);

    EGRID_HD static inline float  eToE_static_log  (int e, int nbins_, float Emin_, float decades_);
    EGRID_HD static inline float  eToEc_static_log (int e, int nbins_, float Emin_, float decades_);
    EGRID_HD static inline float  eToEcg_static_log(int e, int nbins_, float Emin_, float decades_);
    EGRID_HD static inline float  eToEca_static_log(int e, int nbins_, float Emin_, float decades_);
    EGRID_HD static inline float  eTodE_static_log (int e, int nbins_, float Emin_, float decades_);

    // --- Static helpers (linear) ---
    EGRID_HD static inline double eToE_static_linear  (int e, int nbins_, double Emin_, double Emax_);
    EGRID_HD static inline double eToEc_static_linear (int e, int nbins_, double Emin_, double Emax_);
    EGRID_HD static inline double eToEcg_static_linear(int e, int nbins_, double Emin_, double Emax_);
    EGRID_HD static inline double eToEca_static_linear(int e, int nbins_, double Emin_, double Emax_);
    EGRID_HD static inline double eTodE_static_linear (int e, int nbins_, double Emin_, double Emax_);

    EGRID_HD static inline float  eToE_static_linear  (int e, int nbins_, float Emin_, float Emax_);
    EGRID_HD static inline float  eToEc_static_linear (int e, int nbins_, float Emin_, float Emax_);
    EGRID_HD static inline float  eToEcg_static_linear(int e, int nbins_, float Emin_, float Emax_);
    EGRID_HD static inline float  eToEca_static_linear(int e, int nbins_, float Emin_, float Emax_);
    EGRID_HD static inline float  eTodE_static_linear (int e, int nbins_, float Emin_, float Emax_);

    // --- Typed dispatch (decades_or_Emax means: decades for log, Emax for linear) ---
    EGRID_HD static inline float  eToE_static_typed  (int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t);
    EGRID_HD static inline float  eToEc_static_typed (int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t);
    EGRID_HD static inline float  eToEcg_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t);
    EGRID_HD static inline float  eToEca_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t);
    EGRID_HD static inline float  eTodE_static_typed (int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t);

    EGRID_HD static inline double eToE_static_typed  (int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t);
    EGRID_HD static inline double eToEc_static_typed (int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t);
    EGRID_HD static inline double eToEcg_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t);
    EGRID_HD static inline double eToEca_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t);
    EGRID_HD static inline double eTodE_static_typed (int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t);


};

// Logarithmic grid concrete helpers
struct LogEgrid {
    // Edge energies: e=0..nbins
    EGRID_HD static inline double eToE_impl(int e, int nbins, double Emin, double decades) {
        // edges are uniformly spaced in log10 between Emin and Emin*10^decades
        double t = static_cast<double>(e) / static_cast<double>(nbins);
        return Emin * std::pow(10.0, t * decades);
    }

    // Bin centre: geometric mean of edges e and e+1
    EGRID_HD static inline double eToEc_impl(int e, int nbins, double Emin, double decades) {
        double a = eToE_impl(e,     nbins, Emin, decades);
        double b = eToE_impl(e + 1, nbins, Emin, decades);
        return std::sqrt(a * b);
    }

    EGRID_HD static inline double eToEcg_impl(int e, int nbins, double Emin, double decades) {
        return eToEc_impl(e, nbins, Emin, decades);
    }

    // Bin centre: arithmetic mean of edges e and e+1
    EGRID_HD static inline double eToEca_impl(int e, int nbins, double Emin, double decades) {
        double a = eToE_impl(e,     nbins, Emin, decades);
        double b = eToE_impl(e + 1, nbins, Emin, decades);
        return 0.5 * (a + b);
    }

    EGRID_HD static inline double eTodE_impl(int e, int nbins, double Emin, double decades) {
        double a = eToE_impl(e,     nbins, Emin, decades);
        double b = eToE_impl(e + 1, nbins, Emin, decades);
        return b - a;
    }

    // Map energy to bin index using floor (bin edges)
    EGRID_HD static inline int EToIndex_impl(double E, int nbins, double Emin, double decades) {
        if (E <= Emin) return 0;
        double Emax = Emin * std::pow(10.0, decades);
        if (E >= Emax) return nbins - 1;

        double logE = std::log10(E / Emin);         // in [0, decades)
        double x = (logE / decades) * nbins;        // in [0, nbins)
        int idx = static_cast<int>(x);              // floor
        if (idx < 0) idx = 0;
        if (idx >= nbins) idx = nbins - 1;
        return idx;
    }

    EGRID_HD static inline int EToIndex_static(double E, int nbins_, double Emin_, double decades_) {
        return EToIndex_impl(E, nbins_, Emin_, decades_);
    }

    EGRID_HD static inline int EToIndex_static(float E, int nbins_, float Emin_, float decades_) {
        if (E <= Emin_) return 0;
        float Emax_ = Emin_ * std::pow(10.0f, decades_);
        if (E >= Emax_) return nbins_ - 1;

        float logE = std::log10(E / Emin_);
        float x = (logE / decades_) * static_cast<float>(nbins_);
        int idx = static_cast<int>(x); // floor
        if (idx < 0) idx = 0;
        if (idx >= nbins_) idx = nbins_ - 1;
        return idx;
    }

        EGRID_HD static inline float eToE_impl(int e, int nbins, float Emin, float decades) {
        float t = static_cast<float>(e) / static_cast<float>(nbins);
        return Emin * powf(10.0f, t * decades);
    }

    EGRID_HD static inline float eToEc_impl(int e, int nbins, float Emin, float decades) {
        float a = eToE_impl(e,     nbins, Emin, decades);
        float b = eToE_impl(e + 1, nbins, Emin, decades);
        return sqrtf(a * b);
    }

    EGRID_HD static inline float eToEcg_impl(int e, int nbins, float Emin, float decades) {
        return eToEc_impl(e, nbins, Emin, decades);
    }

    EGRID_HD static inline float eToEca_impl(int e, int nbins, float Emin, float decades) {
        float a = eToE_impl(e,     nbins, Emin, decades);
        float b = eToE_impl(e + 1, nbins, Emin, decades);
        return 0.5f * (a + b);
    }

    EGRID_HD static inline float eTodE_impl(int e, int nbins, float Emin, float decades) {
        float a = eToE_impl(e,     nbins, Emin, decades);
        float b = eToE_impl(e + 1, nbins, Emin, decades);
        return b - a;
    }

};

// Linear grid concrete helpers
struct LinEgrid {
    // Edge energies: e=0..nbins
    EGRID_HD static inline double eToE_impl(int e, int nbins, double Emin, double Emax) {
        double t = static_cast<double>(e) / static_cast<double>(nbins);
        return Emin + t * (Emax - Emin);
    }

    // Bin centre: arithmetic mean of edges e and e+1
    EGRID_HD static inline double eToEc_impl(int e, int nbins, double Emin, double Emax) {
        double a = eToE_impl(e,     nbins, Emin, Emax);
        double b = eToE_impl(e + 1, nbins, Emin, Emax);
        return 0.5 * (a + b);
    }

    EGRID_HD static inline double eToEca_impl(int e, int nbins, double Emin, double Emax) {
        return eToEc_impl(e, nbins, Emin, Emax);
    }

    // Bin centre: geometric mean of edges e and e+1
    EGRID_HD static inline double eToEcg_impl(int e, int nbins, double Emin, double Emax) {
        double a = eToE_impl(e,     nbins, Emin, Emax);
        double b = eToE_impl(e + 1, nbins, Emin, Emax);
        return std::sqrt(a * b);
    }

    EGRID_HD static inline double eTodE_impl(int e, int nbins, double Emin, double Emax) {
        double a = eToE_impl(e,     nbins, Emin, Emax);
        double b = eToE_impl(e + 1, nbins, Emin, Emax);
        return b - a;
    }

    // Map energy to bin index using floor (bin edges)
    EGRID_HD static inline int EToIndex_impl(double E, int nbins, double Emin, double Emax) {
        if (E <= Emin) return 0;
        if (E >= Emax) return nbins - 1;

        double t = (E - Emin) / (Emax - Emin);   // in [0,1)
        double x = t * nbins;                    // in [0, nbins)
        int idx = static_cast<int>(x);           // floor
        if (idx < 0) idx = 0;
        if (idx >= nbins) idx = nbins - 1;
        return idx;
    }

    EGRID_HD static inline int EToIndex_static(double E, int nbins_, double Emin_, double Emax_) {
        return EToIndex_impl(E, nbins_, Emin_, Emax_);
    }

    EGRID_HD static inline int EToIndex_static(float E, int nbins_, float Emin_, float Emax_) {
        if (E <= Emin_) return 0;
        if (E >= Emax_) return nbins_ - 1;

        float t = (E - Emin_) / (Emax_ - Emin_);     // in [0,1)
        float x = t * static_cast<float>(nbins_);    // in [0, nbins)
        int idx = static_cast<int>(x);               // floor
        if (idx < 0) idx = 0;
        if (idx >= nbins_) idx = nbins_ - 1;
        return idx;
    }

        EGRID_HD static inline float eToE_impl(int e, int nbins, float Emin, float Emax) {
        float t = static_cast<float>(e) / static_cast<float>(nbins);
        return Emin + t * (Emax - Emin);
    }

    EGRID_HD static inline float eToEc_impl(int e, int nbins, float Emin, float Emax) {
        float a = eToE_impl(e,     nbins, Emin, Emax);
        float b = eToE_impl(e + 1, nbins, Emin, Emax);
        return 0.5f * (a + b);
    }

    EGRID_HD static inline float eToEca_impl(int e, int nbins, float Emin, float Emax) {
        return eToEc_impl(e, nbins, Emin, Emax);
    }

    EGRID_HD static inline float eToEcg_impl(int e, int nbins, float Emin, float Emax) {
        float a = eToE_impl(e,     nbins, Emin, Emax);
        float b = eToE_impl(e + 1, nbins, Emin, Emax);
        return sqrtf(a * b);
    }

    EGRID_HD static inline float eTodE_impl(int e, int nbins, float Emin, float Emax) {
        float a = eToE_impl(e,     nbins, Emin, Emax);
        float b = eToE_impl(e + 1, nbins, Emin, Emax);
        return b - a;
    }


};

// Provide the Egrid instance method definitions now that concrete helpers exist.
inline double Egrid::eToE(int e) const {
    // Clamp edge index to [0, nbins]
    if (e < 0) e = 0;
    if (e > nbins) e = nbins;

    if (type == EgridType::Logarithmic) return LogEgrid::eToE_impl(e, nbins, Emin, decades);
    return LinEgrid::eToE_impl(e, nbins, Emin, Emax);
}

inline double Egrid::eToEcg(int e) const {
    // Clamp bin index to [0, nbins-1]
    if (e < 0) e = 0;
    if (e >= nbins) e = nbins - 1;

    if (type == EgridType::Logarithmic) return LogEgrid::eToEcg_impl(e, nbins, Emin, decades);
    return LinEgrid::eToEcg_impl(e, nbins, Emin, Emax);
}

inline double Egrid::eToEca(int e) const {
    // Clamp bin index to [0, nbins-1]
    if (e < 0) e = 0;
    if (e >= nbins) e = nbins - 1;

    if (type == EgridType::Logarithmic) return LogEgrid::eToEca_impl(e, nbins, Emin, decades);
    return LinEgrid::eToEca_impl(e, nbins, Emin, Emax);
}

inline double Egrid::eToEc(int e) const { return eToEcg(e); }

inline double Egrid::eTodE(int e) const {
    // Clamp bin index to [0, nbins-1]
    if (e < 0) e = 0;
    if (e >= nbins) e = nbins - 1;

    if (type == EgridType::Logarithmic) return LogEgrid::eTodE_impl(e, nbins, Emin, decades);
    return LinEgrid::eTodE_impl(e, nbins, Emin, Emax);
}

inline std::vector<double> Egrid::getEdges() const {
    std::vector<double> edges(static_cast<std::size_t>(nbins) + 1);
    for (int e = 0; e <= nbins; ++e) {
        edges[static_cast<std::size_t>(e)] = eToE(e);
    }
    return edges;
}

inline int Egrid::EToIndex(double E) const {
    if (type == EgridType::Logarithmic) return LogEgrid::EToIndex_impl(E, nbins, Emin, decades);
    return LinEgrid::EToIndex_impl(E, nbins, Emin, Emax);
}

// Backwards-compatible static alias uses logarithmic behaviour
EGRID_HD inline int Egrid::EToIndex_static(double E, int nbins_, double Emin_, double decades_) {
    return LogEgrid::EToIndex_static(E, nbins_, Emin_, decades_);
}

EGRID_HD inline int Egrid::EToIndex_static(float E, int nbins_, float Emin_, float decades_) {
    return LogEgrid::EToIndex_static(E, nbins_, Emin_, decades_);
}


EGRID_HD inline int Egrid::EToIndex_static_linear(double E, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::EToIndex_static(E, nbins_, Emin_, Emax_);
}

EGRID_HD inline int Egrid::EToIndex_static_linear(float E, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::EToIndex_static(E, nbins_, Emin_, Emax_);
}

EGRID_HD inline double Egrid::eToE_static_log(int e, int nbins_, double Emin_, double decades_) {
    return LogEgrid::eToE_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline double Egrid::eToEc_static_log(int e, int nbins_, double Emin_, double decades_) {
    return LogEgrid::eToEc_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline double Egrid::eToEcg_static_log(int e, int nbins_, double Emin_, double decades_) {
    return LogEgrid::eToEcg_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline double Egrid::eToEca_static_log(int e, int nbins_, double Emin_, double decades_) {
    return LogEgrid::eToEca_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline double Egrid::eTodE_static_log(int e, int nbins_, double Emin_, double decades_) {
    return LogEgrid::eTodE_impl(e, nbins_, Emin_, decades_);
}

EGRID_HD inline float Egrid::eToE_static_log(int e, int nbins_, float Emin_, float decades_) {
    return LogEgrid::eToE_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline float Egrid::eToEc_static_log(int e, int nbins_, float Emin_, float decades_) {
    return LogEgrid::eToEc_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline float Egrid::eToEcg_static_log(int e, int nbins_, float Emin_, float decades_) {
    return LogEgrid::eToEcg_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline float Egrid::eToEca_static_log(int e, int nbins_, float Emin_, float decades_) {
    return LogEgrid::eToEca_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline float Egrid::eTodE_static_log(int e, int nbins_, float Emin_, float decades_) {
    return LogEgrid::eTodE_impl(e, nbins_, Emin_, decades_);
}
EGRID_HD inline double Egrid::eToE_static_linear(int e, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::eToE_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline double Egrid::eToEc_static_linear(int e, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::eToEc_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline double Egrid::eToEcg_static_linear(int e, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::eToEcg_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline double Egrid::eToEca_static_linear(int e, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::eToEca_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline double Egrid::eTodE_static_linear(int e, int nbins_, double Emin_, double Emax_) {
    return LinEgrid::eTodE_impl(e, nbins_, Emin_, Emax_);
}

EGRID_HD inline float Egrid::eToE_static_linear(int e, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::eToE_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline float Egrid::eToEc_static_linear(int e, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::eToEc_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline float Egrid::eToEcg_static_linear(int e, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::eToEcg_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline float Egrid::eToEca_static_linear(int e, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::eToEca_impl(e, nbins_, Emin_, Emax_);
}
EGRID_HD inline float Egrid::eTodE_static_linear(int e, int nbins_, float Emin_, float Emax_) {
    return LinEgrid::eTodE_impl(e, nbins_, Emin_, Emax_);
}


// Typed dispatch: decades_or_Emax_ means:
//   - if t==Logarithmic: parameter is decades
//   - if t==Linear:      parameter is Emax
EGRID_HD inline int Egrid::EToIndex_static_typed(double E, int nbins_, double Emin_, double decades_or_Emax_,
                                                EgridType t) {
    if (t == EgridType::Logarithmic) {
        return LogEgrid::EToIndex_static(E, nbins_, Emin_, decades_or_Emax_);
    } else {
        return LinEgrid::EToIndex_static(E, nbins_, Emin_, decades_or_Emax_);
    }
}

EGRID_HD inline int Egrid::EToIndex_static_typed(float E, int nbins_, float Emin_, float decades_or_Emax_,
                                                EgridType t) {
    if (t == EgridType::Logarithmic) {
        return LogEgrid::EToIndex_static(E, nbins_, Emin_, decades_or_Emax_);
    } else {
        return LinEgrid::EToIndex_static(E, nbins_, Emin_, decades_or_Emax_);
    }
}

EGRID_HD inline double Egrid::eToE_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToE_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToE_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline double Egrid::eToEc_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEc_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEc_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline double Egrid::eToEcg_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEcg_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEcg_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline double Egrid::eToEca_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEca_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEca_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline double Egrid::eTodE_static_typed(int e, int nbins_, double Emin_, double decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eTodE_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eTodE_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}

EGRID_HD inline float Egrid::eToE_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToE_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToE_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline float Egrid::eToEc_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEc_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEc_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline float Egrid::eToEcg_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEcg_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEcg_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline float Egrid::eToEca_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eToEca_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eToEca_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}
EGRID_HD inline float Egrid::eTodE_static_typed(int e, int nbins_, float Emin_, float decades_or_Emax_, EgridType t) {
    return (t == EgridType::Logarithmic)
        ? eTodE_static_log(e, nbins_, Emin_, decades_or_Emax_)
        : eTodE_static_linear(e, nbins_, Emin_, decades_or_Emax_);
}

#endif // EGRID_HPP
