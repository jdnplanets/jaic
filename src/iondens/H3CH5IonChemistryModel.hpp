/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
 The H3CH5IonChemistryModel implements a steady-state, no-transport chemistry system involving H3+ and CH5+ ions, 
with production and loss processes defined by ionization rates and reaction rates.
 */

 // H3CH5IonChemistryModel.hpp

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "IonChemistryModel.hpp"
#include "IonChemistryRates.hpp"

// Steady-state, no-transport two-ion chemistry:
//   H3+ produced at total rate q (via rapid H2+ + H2 -> H3+)
//   H3+ destroyed by dissociative recombination and reaction with methane producing CH5+
//   CH5+ destroyed by dissociative recombination
//
// Equations (densities in cm^-3, rates in cm^3 s^-1, q in cm^-3 s^-1):
//   0 = q - alpha_H3 * nH3 * (nH3 + nCH5) - k * nH3 * nH4
//   0 = k * nH3 * nCH4 - alpha_CH5 * nCH5 * (nH3 + nCH5)
//
// We solve a single cubic for S = nH3 + nCH5 via robust bisection, then back-substitute.
struct H3CH5IonChemistryModel final : IIonChemistryModel {
    static double fS(double S, double q, double a1, double a2, double km) {
        // f(S) = (a2*a1) S^3 + (a2*km) S^2 - (q*a2) S - (q*km)
        return (a2 * a1) * S * S * S + (a2 * km) * S * S - (q * a2) * S - (q * km);
    }

    static double solve_positive_root_bisection(double q, double a1, double a2, double km) {
        // For q>0 and km>0:
        //   f(0) = -q*km < 0
        //   f(S) -> +inf as S->+inf
        // so there is at least one positive root; bisection is safe if we bracket it.
        double lo = 0.0;
        double hi = 1.0;

        // Grow hi until f(hi) > 0
        for (int i = 0; i < 200; ++i) {
            if (fS(hi, q, a1, a2, km) > 0.0) break;
            hi *= 2.0;
            if (hi > 1e30) break;
        }

        // If we failed to bracket, return the last hi.
        if (fS(hi, q, a1, a2, km) <= 0.0) return hi;

        for (int it = 0; it < 120; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fm  = fS(mid, q, a1, a2, km);
            if (fm > 0.0) hi = mid;
            else          lo = mid;
        }
        return 0.5 * (lo + hi);
    }

    void compute(const Params&,
                 World1D& world,
                 const std::vector<double>& Qtot) const override
    {
        for (int z = 0; z < world.nz; ++z) {
            const double q = std::max(0.0, Qtot[z]);
            const double nCH4 = (z < static_cast<int>(world.nCH4.size())) ? std::max(0.0f, world.nCH4[z]) : 0.0;
            const double T = (z < static_cast<int>(world.T.size())) ? world.T[z] : 0.0;

            const double a1 = IonChemistryRates::alpha_H3p(T);
            const double a2 = IonChemistryRates::alpha_CH5p(T);
            const double k  = IonChemistryRates::k_H3p_CH4;

            const double km = k * nCH4;

            double nH3 = 0.0;
            double nCH5 = 0.0;

            if (q <= 0.0) {
                nH3 = 0.0;
                nCH5 = 0.0;
            } else if (km <= 0.0 || a2 <= 0.0) {
                // No methane conversion or invalid a2 -> reduces to pure H3+ balance.
                nH3 = (a1 > 0.0) ? std::sqrt(q / a1) : 0.0;
                nCH5 = 0.0;
            } else if (a1 <= 0.0) {
                // Invalid a1; keep zeros rather than NaN.
                nH3 = 0.0;
                nCH5 = 0.0;
            } else {
                // Solve for S = nH3 + nCH5.
                const double S = solve_positive_root_bisection(q, a1, a2, km);

                // Back-substitution: nH3 = (a2 S^2) / (km + a2 S)
                nH3 = (a2 * S * S) / (km + a2 * S);
                nCH5 = S - nH3;

                if (nH3 < 0.0) nH3 = 0.0;
                if (nCH5 < 0.0) nCH5 = 0.0;
            }

            world.nH3p[z]   = static_cast<float>(nH3);
            world.nCH5p[z]  = static_cast<float>(nCH5);
            world.nC3Hnp[z] = 0.0f; // not included in this two-ion closure
            world.ne[z]     = world.nH3p[z] + world.nCH5p[z] + world.nC3Hnp[z];
        }
    }
};
