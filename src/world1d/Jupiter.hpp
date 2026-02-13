/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the Jupiter class, which inherits from World1D and populates the atmospheric state 
(nH2, T, hydrocarbons) based on a combination of binary files (for nH2) and CSV files (for hydrocarbon mixing ratios). 
The temperature profile is a simple piecewise linear ramp defined by parameters in the Config struct. 
*/

 // Jupiter.hpp

#pragma once
#include "World1D.hpp"
#include "WorldOverrides.hpp"
#include "FileIO.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>     // getenv
#include <stdexcept>
#include <utility>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

struct Jupiter final : public World1D {

    struct HydrocarbonMixingRatios {
        double mix = 0.0;  // mixing ratio
        double z_km = 0.0; // altitude in km at which this mixing ratio applies
    };

    struct IonMixingRatios {
        double z;       // altitude in km
        double CH5p;    // CH5+ mixing ratio
        double C3Hnp;   // C3Hn+ mixing ratio
        double H3p;     // H3+ mixing ratio
    };


    struct Config {
        // Domain
        double Z0_m = 0.0;          // bottom of domain (m)
        double Z1_m = 3000e3;       // top of domain (m)
        int    nz   = 300;          // number of altitude bins

        // Magnetic field
        double B_T  = 0.0015;       // Tesla

        // Temperature profile 
        double T0_K     = 200.0;    // temperature at zTmin_km (K)
        double T1_K     = 1200.0;   // temperature at zTmax_km (K)
        double zTmin_km = 200.0;    // altitude where T = T0_K (km)
        double zTmax_km = 1000.0;   // altitude where T = T1_K (km)

        // Hydrocarbon mixing ratio files (CSV: mix, pressure[mbar])
        std::string mixCH4_file  = "mixing_ratio_ch4.csv";
        std::string mixC2H2_file = "mixing_ratio_c2h2.csv";
        std::string mixC2H4_file = "mixing_ratio_c2h4.csv";
        std::string mixC2H6_file = "mixing_ratio_c2h6.csv";

        // H2 density binary file
        std::string nH2_file = "nH2vz.bin";

        // Data directory
        std::string data_dir;

        // Pressure->altitude conversion parameters 
        double H_km  = 27.0;    // scale height (km)
        double p0_mbar = 1e3;   // reference pressure (mbar)
    };

    static Config Defaults() {
        Config c;
        c.data_dir = datadir;
        return c;
    }

    // Default Jupiter world (uses Defaults()).
    Jupiter()
        : Jupiter(Defaults()) {}

    // Construct with custom config if desired.
    explicit Jupiter(const Config& cfg_)
        : World1D(static_cast<float>(cfg_.Z0_m),
                  static_cast<float>(cfg_.Z1_m),
                  cfg_.nz),
          cfg(cfg_)
    {
        B = cfg.B_T;
    }

    // Allow the generic override bag to adjust B, etc., before init().
    void apply_overrides(const WorldOverrides& o) override {
        World1D::apply_overrides(o);
        // Keep cfg consistent with base
        cfg.B_T = B;
    }


protected:
    void populate() override {
        // nH2 from binary file
        read_nH2_binary();

        // Temperature profile (simple piecewise linear ramp)
        fill_temperature_profile();

        // Hydrocarbon densities from mixing ratio CSVs
        //    density(z) = nH2(z) * mix(z)
        nCH4  = compute_hydrocarbon_density(cfg.mixCH4_file);
        nC2H2 = compute_hydrocarbon_density(cfg.mixC2H2_file);
        nC2H4 = compute_hydrocarbon_density(cfg.mixC2H4_file);
        nC2H6 = compute_hydrocarbon_density(cfg.mixC2H6_file);

    }

    void validate() const override {
        World1D::validate();
        if (cfg.nz != nz) throw std::runtime_error("Jupiter: cfg.nz inconsistent with base nz");
        if (cfg.Z1_m <= cfg.Z0_m) throw std::runtime_error("Jupiter: Z1_m must be > Z0_m");
        if (cfg.data_dir.empty()) throw std::runtime_error("Jupiter: data_dir is empty");
        if (cfg.H_km <= 0.0) throw std::runtime_error("Jupiter: H_km must be > 0");
        if (cfg.p0_mbar <= 0.0) throw std::runtime_error("Jupiter: p0_mbar must be > 0");
    }

private:
    Config cfg;



    // ---------------------------
    // nH2 binary reader
    // ---------------------------
    void read_nH2_binary() {
        const std::string fullpath = path_join(cfg.data_dir, cfg.nH2_file);
        utils::io::readFloatBinaryfileData(fullpath, nH2.data(), nz, 1);
    }

    // ---------------------------
    // Temperature profile
    // ---------------------------
    void fill_temperature_profile() {
        const double T0 = cfg.T0_K;
        const double T1 = cfg.T1_K;
        const double zmin = cfg.zTmin_km;
        const double zmax = cfg.zTmax_km;

        const double m = (zmax != zmin) ? (T1 - T0) / (zmax - zmin) : 0.0;

        for (int i = 0; i < nz; ++i) {
            const double z_km = static_cast<double>(Z[i]) / 1e3;

            if (z_km < zmin) {
                T[i] = T0;
            } else if (z_km < zmax) {
                T[i] = m * (z_km - zmin) + T0;
            } else {
                T[i] = T1;
            }
        }
    }

    // ---------------------------
    // Mixing ratio CSV helpers
    // CSV columns: mix, pressure(mbar)
    // Converts pressure -> altitude via z = H ln(p0/p)
    // ---------------------------
    std::vector<HydrocarbonMixingRatios>
    read_hydrocarbon_mixing_ratios(const std::string& filename) const {
        const std::string fullpath = path_join(cfg.data_dir, filename);
        std::ifstream file(fullpath);
        if (!file) {
            throw std::runtime_error("Jupiter: failed to open mixing ratio file: " + fullpath);
        }

        std::vector<HydrocarbonMixingRatios> data;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            std::string a, b;

            if (!std::getline(iss, a, ',')) continue;
            if (!std::getline(iss, b, ',')) continue;

            const double mix = std::stod(a);
            const double p_mbar = std::stod(b);

            if (p_mbar <= 0.0) continue; // skip invalid

            HydrocarbonMixingRatios pt;
            pt.mix = mix;
            pt.z_km = cfg.H_km * std::log(cfg.p0_mbar / p_mbar);

            data.push_back(pt);
        }

        if (data.empty()) {
            throw std::runtime_error("Jupiter: mixing ratio file had no usable rows: " + fullpath);
        }

        // Ensure sorted by altitude (in case file isn't)
        std::sort(data.begin(), data.end(),
                  [](const HydrocarbonMixingRatios& x, const HydrocarbonMixingRatios& y) {
                      return x.z_km < y.z_km;
                  });

        return data;
    }

    static double get_mix_at_z_km(double z_km,
                                 const std::vector<HydrocarbonMixingRatios>& table)
    {
        if (z_km <= table.front().z_km) return table.front().mix;
        if (z_km >= table.back().z_km)  return table.back().mix; // last value tiny anyway

        // linear interpolation in z
        for (std::size_t i = 1; i < table.size(); ++i) {
            if (z_km < table[i].z_km) {
                const double z1 = table[i - 1].z_km;
                const double z2 = table[i].z_km;
                const double m1 = table[i - 1].mix;
                const double m2 = table[i].mix;
                return m1 + (m2 - m1) * (z_km - z1) / (z2 - z1);
            }
        }
        return table.back().mix;
    }

    std::vector<double> compute_hydrocarbon_density(const std::string& mix_file) const {
        const auto mix_table = read_hydrocarbon_mixing_ratios(mix_file);

        std::vector<double> density(static_cast<std::size_t>(nz), 0.0);
        for (int i = 0; i < nz; ++i) {
            const double z_km = static_cast<double>(Z[i]) / 1e3;
            const double mix = get_mix_at_z_km(z_km, mix_table);
            density[static_cast<std::size_t>(i)] = static_cast<double>(nH2[i]) * mix;
        }
        return density;
    }



};
