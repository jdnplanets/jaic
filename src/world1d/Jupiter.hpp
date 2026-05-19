/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the Jupiter class, which inherits from World1D and populates the atmospheric state
from an ASCII atmosphere file containing pressure, altitude, temperature, and neutral number densities.
*/

// Jupiter.hpp

#pragma once

#include "World1D.hpp"
#include "WorldOverrides.hpp"
#include "FileIO.hpp"
#include "ArrayOps.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

struct Jupiter final : public World1D {

    struct Config {

        // Atmosphere definition file
        std::string data_dir;
        std::string atm_file = "atmosphere.dat";

        // Domain
        double P0 = 1e5;   // bottom of domain, Pa
        double P1 = 1e-7;  // top of domain, Pa
        int    nP = 500;   // number of pressure levels

        // Magnetic field
        double B_T = 0.001;
        double Bdipang_deg = 90.0;
    };

    static Config Defaults() {
        Config c;
        c.data_dir = datadir;
        return c;
    }

    Jupiter()
        : Jupiter(Defaults()) {}

    explicit Jupiter(const Config& cfg_)
        : World1D(static_cast<float>(cfg_.P0),
                  static_cast<float>(cfg_.P1),
                  cfg_.nP),
          cfg(cfg_)
    {
        B = cfg.B_T;
        Bdipang = cfg.Bdipang_deg;
    }

    void apply_overrides(const WorldOverrides& o) override {
        World1D::apply_overrides(o);

        cfg.B_T = B;
        cfg.Bdipang_deg = Bdipang;
    }

protected:

    void populate() override {

        const std::string filepath = path_join(cfg.data_dir, cfg.atm_file);

        AtmosphereFile atm = read_atmosphere_file(filepath);

        validate_atmosphere_file(atm, filepath);

        /*
         * If later we want cfg.nP to differ from the atmosphere file length,
         * then this is where we would interpolate atm.P, atm.Z, atm.T, etc.
         * onto the model pressure grid.
         */
        if (static_cast<int>(atm.P.size()) != nP) {
            throw std::runtime_error(
                "Jupiter: atmosphere file has " +
                std::to_string(atm.P.size()) +
                " rows, but World1D expects nP = " +
                std::to_string(nP) +
                ". Add pressure-grid interpolation here."
            );
        }

        for (int i = 0; i < nP; ++i) {
            P[i]     = static_cast<float>(atm.P[i]);      // Pa
            Z[i]     = static_cast<float>(atm.Z[i]);      // m
            T[i]     = static_cast<float>(atm.T[i]);      // K

            nH2[i]   = static_cast<float>(atm.nH2[i]);    // cm^-3
            nHe[i]   = static_cast<float>(atm.nHe[i]);
            nH[i]    = static_cast<float>(atm.nH[i]);
            nCH4[i]  = static_cast<float>(atm.nCH4[i]);
            nC2H2[i] = static_cast<float>(atm.nC2H2[i]);
            nC2H4[i] = static_cast<float>(atm.nC2H4[i]);
            nC2H6[i] = static_cast<float>(atm.nC2H6[i]);
        }

        compute_grid_spacings_from_Z();
    }

    void validate() const override {
        World1D::validate();

        if (cfg.nP != nP) {
            throw std::runtime_error("Jupiter: cfg.nP inconsistent with base nP");
        }
        if (!std::isfinite(cfg.P0) || cfg.P0 <= 0.0) {
            throw std::runtime_error("Jupiter: cfg.P0 must be finite and > 0 Pa");
        }
        if (!std::isfinite(cfg.P1) || cfg.P1 <= 0.0) {
            throw std::runtime_error("Jupiter: cfg.P1 must be finite and > 0 Pa");
        }
        if (cfg.P0 <= cfg.P1) {
            throw std::runtime_error(
                "Jupiter: expected cfg.P0 > cfg.P1 for a pressure grid running upward"
            );
        }
        if (cfg.nP < 2) {
            throw std::runtime_error("Jupiter: cfg.nP must be at least 2");
        }
        if (cfg.data_dir.empty()) {
            throw std::runtime_error("Jupiter: cfg.data_dir is empty");
        }
        if (cfg.atm_file.empty()) {
            throw std::runtime_error("Jupiter: cfg.atm_file is empty");
        }
        if (!std::isfinite(cfg.B_T) || cfg.B_T <= 0.0) {
            throw std::runtime_error("Jupiter: cfg.B_T must be finite and > 0");
        }
        if (!std::isfinite(cfg.Bdipang_deg)) {
            throw std::runtime_error("Jupiter: cfg.Bdipang_deg must be finite");
        }
        if (cfg.Bdipang_deg < 0.0 || cfg.Bdipang_deg > 90.0) {
            throw std::runtime_error("Jupiter: cfg.Bdipang_deg must be in [0, 90] degrees");
        }
        validate_populated_grid();
    }

private:

    struct AtmosphereFile {
        std::vector<double> P;      // Pa, converted from bar in file
        std::vector<double> Z;      // m, converted from km in file
        std::vector<double> T;      // K

        std::vector<double> nH2;    // cm^-3
        std::vector<double> nHe;    // cm^-3
        std::vector<double> nH;     // cm^-3
        std::vector<double> nCH4;   // cm^-3
        std::vector<double> nC2H2;  // cm^-3
        std::vector<double> nC2H4;  // cm^-3
        std::vector<double> nC2H6;  // cm^-3
    };

    Config cfg;

    static std::string trim(const std::string& s) {
        const auto first = s.find_first_not_of(" \t\r\n");

        if (first == std::string::npos) {
            return "";
        }

        const auto last = s.find_last_not_of(" \t\r\n");
        return s.substr(first, last - first + 1);
    }

    static std::vector<double> parse_line(const std::string& line) {
        std::vector<double> values;

        std::stringstream ss(line);
        double x;

        while (ss >> x) {
            values.push_back(x);
        }

        if (!ss.eof()) {
            throw std::runtime_error(
                "failed to parse whitespace-separated numeric line: " + line
            );
        }

        return values;
    }

    static AtmosphereFile read_atmosphere_file(const std::string& filepath) {
        std::ifstream file(filepath);

        if (!file) {
            throw std::runtime_error(
                "Jupiter: could not open atmosphere file: " + filepath
            );
        }

        AtmosphereFile atm;

        std::string line;
        std::size_t line_number = 0;

        while (std::getline(file, line)) {
            ++line_number;

            line = trim(line);

            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::vector<double> v;

            try {
                v = parse_line(line);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "Jupiter: failed to parse atmosphere file line " +
                    std::to_string(line_number) + ": " + e.what()
                );
            }

            if (v.size() != 10) {
                throw std::runtime_error(
                    "Jupiter: expected 10 columns in atmosphere file line " +
                    std::to_string(line_number) +
                    ", found " + std::to_string(v.size())
                );
            }

            const double p_bar = v[0];
            const double z_km  = v[1];

            atm.P.push_back(p_bar * 1.0e5);  // bar -> Pa
            atm.Z.push_back(z_km  * 1.0e3);  // km  -> m
            atm.T.push_back(v[2]);

            atm.nH2.push_back(v[3]);
            atm.nHe.push_back(v[4]);
            atm.nH.push_back(v[5]);
            atm.nCH4.push_back(v[6]);
            atm.nC2H2.push_back(v[7]);
            atm.nC2H4.push_back(v[8]);
            atm.nC2H6.push_back(v[9]);
        }

        return atm;
    }

    void validate_atmosphere_file(
        const AtmosphereFile& atm,
        const std::string& filepath
    ) const {
        const std::size_t n = atm.P.size();

        if (n < 2) {
            throw std::runtime_error(
                "Jupiter: atmosphere file must contain at least two data rows: " +
                filepath
            );
        }

        const auto same_size = [n](const std::vector<double>& x) {
            return x.size() == n;
        };

        if (!same_size(atm.Z)     ||
            !same_size(atm.T)     ||
            !same_size(atm.nH2)   ||
            !same_size(atm.nHe)   ||
            !same_size(atm.nH)    ||
            !same_size(atm.nCH4)  ||
            !same_size(atm.nC2H2) ||
            !same_size(atm.nC2H4) ||
            !same_size(atm.nC2H6)) {
            throw std::runtime_error(
                "Jupiter: inconsistent atmosphere column lengths in " + filepath
            );
        }

        for (std::size_t i = 0; i < n; ++i) {
            check_finite_positive(atm.P[i], "P", i);
            check_finite(atm.Z[i], "Z", i);
            check_finite_positive(atm.T[i], "T", i);

            check_finite_nonnegative(atm.nH2[i],   "nH2", i);
            check_finite_nonnegative(atm.nHe[i],   "nHe", i);
            check_finite_nonnegative(atm.nH[i],    "nH", i);
            check_finite_nonnegative(atm.nCH4[i],  "nCH4", i);
            check_finite_nonnegative(atm.nC2H2[i], "nC2H2", i);
            check_finite_nonnegative(atm.nC2H4[i], "nC2H4", i);
            check_finite_nonnegative(atm.nC2H6[i], "nC2H6", i);
        }

        /*
         * The atmosphere file should be ordered from deep atmosphere to upper atmosphere:
         *
         *   P decreases with row number
         *   Z increases with row number
         */
        for (std::size_t i = 1; i < n; ++i) {
            if (!(atm.P[i] < atm.P[i - 1])) {
                throw std::runtime_error(
                    "Jupiter: atmosphere pressure must be strictly decreasing with row number; "
                    "problem near row " + std::to_string(i)
                );
            }

            if (!(atm.Z[i] > atm.Z[i - 1])) {
                throw std::runtime_error(
                    "Jupiter: atmosphere altitude must be strictly increasing with row number; "
                    "problem near row " + std::to_string(i)
                );
            }
        }

        const double file_P_bottom = atm.P.front();
        const double file_P_top    = atm.P.back();

        if (cfg.P0 > file_P_bottom) {
            throw std::runtime_error(
                "Jupiter: requested cfg.P0 = " + std::to_string(cfg.P0) +
                " Pa is deeper than the atmosphere file bottom pressure = " +
                std::to_string(file_P_bottom) + " Pa"
            );
        }

        if (cfg.P1 < file_P_top) {
            throw std::runtime_error(
                "Jupiter: requested cfg.P1 = " + std::to_string(cfg.P1) +
                " Pa is above the atmosphere file top pressure = " +
                std::to_string(file_P_top) + " Pa"
            );
        }
    }

    void validate_populated_grid() const {
        for (int i = 0; i < nP; ++i) {
            check_finite_positive(P[i], "P", static_cast<std::size_t>(i));
            check_finite(Z[i], "Z", static_cast<std::size_t>(i));
            check_finite_positive(T[i], "T", static_cast<std::size_t>(i));

            check_finite_nonnegative(nH2[i],   "nH2", static_cast<std::size_t>(i));
            check_finite_nonnegative(nHe[i],   "nHe", static_cast<std::size_t>(i));
            check_finite_nonnegative(nH[i],    "nH", static_cast<std::size_t>(i));
            check_finite_nonnegative(nCH4[i],  "nCH4", static_cast<std::size_t>(i));
            check_finite_nonnegative(nC2H2[i], "nC2H2", static_cast<std::size_t>(i));
            check_finite_nonnegative(nC2H4[i], "nC2H4", static_cast<std::size_t>(i));
            check_finite_nonnegative(nC2H6[i], "nC2H6", static_cast<std::size_t>(i));

            check_finite_nonnegative(dz[i],   "dz", static_cast<std::size_t>(i));
            check_finite_nonnegative(dzcm[i], "dzcm", static_cast<std::size_t>(i));
        }

        for (int i = 1; i < nP; ++i) {
            if (!(P[i] < P[i - 1])) {
                throw std::runtime_error(
                    "Jupiter: populated pressure grid must be strictly decreasing; "
                    "problem near index " + std::to_string(i)
                );
            }

            if (!(Z[i] > Z[i - 1])) {
                throw std::runtime_error(
                    "Jupiter: populated altitude grid must be strictly increasing; "
                    "problem near index " + std::to_string(i)
                );
            }
        }
    }

    void compute_grid_spacings_from_Z() {
        if (nP < 2) {
            throw std::runtime_error("Jupiter: cannot compute dz with nP < 2");
        }

        for (int i = 0; i < nP; ++i) {
            if (i == 0) {
                dz[i] = Z[1] - Z[0];
            } else if (i == nP - 1) {
                dz[i] = Z[nP - 1] - Z[nP - 2];
            } else {
                dz[i] = 0.5f * (Z[i + 1] - Z[i - 1]);
            }

            dzcm[i] = dz[i] * 100.0f;
        }

        Z0 = Z.front();
        Z1 = Z.back();
        nz = nP;
    }

    template <typename T>
    static void check_finite(T x, const std::string& name, std::size_t i) {
        if (!std::isfinite(static_cast<double>(x))) {
            throw std::runtime_error(
                "Jupiter: non-finite " + name +
                " at index/row " + std::to_string(i)
            );
        }
    }

    template <typename T>
    static void check_finite_positive(T x, const std::string& name, std::size_t i) {
        check_finite(x, name, i);

        if (!(static_cast<double>(x) > 0.0)) {
            throw std::runtime_error(
                "Jupiter: expected " + name +
                " > 0 at index/row " + std::to_string(i)
            );
        }
    }

    template <typename T>
    static void check_finite_nonnegative(T x, const std::string& name, std::size_t i) {
        check_finite(x, name, i);

        if (static_cast<double>(x) < 0.0) {
            throw std::runtime_error(
                "Jupiter: expected " + name +
                " >= 0 at index/row " + std::to_string(i)
            );
        }
    }
};