/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // IonDensity.cpp

#include "IonDensity.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdlib>     // getenv
#include <stdexcept>
#include <cmath>

#include "NullIonChemistryModel.hpp"
#include "FileIO.hpp"

#include "constants.h"
using namespace constants;

static std::string home_dir() {
    const char* h = std::getenv("HOME");
    return h ? std::string(h) : std::string();
}

IonDensity::IonDensity(Params params_,
                       std::shared_ptr<Source> src_,
                       World1D &world_)
    : sp(params_),
      src(std::move(src_)),
      world(world_)
{
    world.nH3p.assign(world.nz, 0.0f);
    world.nCH5p.assign(world.nz, 0.0f);
    world.nC3Hnp.assign(world.nz, 0.0f);
    world.ne.assign(world.nz, 0.0f);

    Qtot.assign(world.nz, 0.0);

    outdir = home_dir() + "/cpp/jaic/out/iondens/" + sp.runid + "/";

    if (!std::filesystem::exists(outdir)) {
        if (!std::filesystem::create_directories(outdir)) {
            throw std::runtime_error("Failed to create output directory: " + outdir);
        }
    }
}


void IonDensity::readIonisationRates(Params& params,
                                     std::shared_ptr<Source> src_in,
                                     float F)
{
    const std::string basepath = "out/precip/" + params.runid + "/";
    const std::string suf = src_in->label(sp);
    const std::string filename = basepath + "qion_" + suf + ".dat";
    std::cout << "Reading ionisation rates from: " << filename << "\n";
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open ionization rates file: " + filename);
    }

    std::string line;
    int z = 0;
    bool isNewFormat = false;

    // Peek at the first non-comment line to determine format
    while (std::getline(infile, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') {
            isNewFormat = true;
            continue;
        }
        // First data line found, rewind and reprocess
        infile.clear();
        infile.seekg(0);
        break;
    }

    if (isNewFormat) {
        std::cout << "Detected new ionization rates file format (with z column)\n";
    } else {
        std::cout << "Detected old ionization rates file format (no z column)\n";
    }

    while (std::getline(infile, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;

        std::istringstream iss(line);

        if (isNewFormat) {
            double P = 0.0;
            double Zkm = 0.0;
            double q = 0.0;
            double q1 = 0.0; // contribution from primaries only. Ignored here.
            if (!(iss >> P >> Zkm >> q >> q1)) continue;
            if (Zkm < world.Z0 || Zkm >= world.Z1) {
                throw std::runtime_error(
                    "Ionisation rates file has out-of-range z_index=" + std::to_string(Zkm) +
                    " (world.Z0=" + std::to_string(world.Z0) + ", world.Z1=" + std::to_string(world.Z1) + ") in " + filename
                );
            }
            Qtot[z] = static_cast<float>(q);
        } else {
            double q = 0.0;
            if (!(iss >> q)) continue;
            Qtot[z] = static_cast<float>(q);
        }
        ++z;
    }

    if (z == 0) {
        throw std::runtime_error("No ionisation rate data rows found in: " + filename);
    }

    for (int z = 0; z < world.nz; ++z) {
        Qtot[z] *= F;
    }
}


void IonDensity::setIonisationRates(std::vector<float>& Qzin, const float F) {
    if (Qzin.size() < static_cast<size_t>(world.nz)) {
        throw std::runtime_error("IonDensity::setIonisationRates: input Qzin too short");
    }
    for (int z = 0; z < world.nz; ++z) {
        Qtot[z] = static_cast<double>(Qzin[z]) * static_cast<double>(F);
        // std::cout << Qtot[z] << "\n";
    }
}

void IonDensity::computeIonDensity() {
    // Ion density is computed by the chemistry model attached to the World.
    // If a world has no model, fall back to the simplest assumption: all H3+.
    if (!world.ionChemistryModel) {
        NullIonChemistryModel fallback;
        fallback.compute(sp, world, Qtot);
        return;
    }

    world.ionChemistryModel->compute(sp, world, Qtot);
}

void IonDensity::writeIonDensityToFile() const {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    const std::string suf = src->label(sp);
    const std::string filename = outdir + "iondens_" + suf + ".dat";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "ion_densities"},
        {"layout", "rows=P_index (0..nP-1)"},
        {"nP", std::to_string(world.nP)},
        {"P_units", "Pa"},
        {"Qtot_units", "cm-3 s-1"},
        {"density_units", "cm-3"},
    };

    const std::vector<std::string> cols = {
        "P [Pa]",
        "Z [km]",
        "Qtot [cm-3_s-1]",
        "nH3p [cm-3]",
        "nCH5p [cm-3]",
        "nC3Hnp [cm-3]",
        "ne [cm-3]"
    };

    const int colw = 14;
    const int precision = 8;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "IonDens model output: ion densities vs altitude",
        meta,
        cols,
        world.nz,
        [&](int z, std::ostream& os, int w) {
            os << std::right
            << std::setw(w) << world.P[z]       << " "
            << std::setw(w) << world.Z[z]/1e3   << " "
            << std::setw(w) << Qtot[z]          << " "
            << std::setw(w) << world.nH3p[z]    << " "
            << std::setw(w) << world.nCH5p[z]   << " "
            << std::setw(w) << world.nC3Hnp[z]  << " "
            << std::setw(w) << world.ne[z];
        },
        colw,
        precision
    );

    if (ok) std::cout << "Ion density written to: " << filename << "\n";
    else    std::cerr << "Failed to write ion density to: " << filename << "\n";
}

double IonDensity::getH3pColumnDensity() const {
    double columnDensity = 0.0;
    for (int z = 0; z < world.nz; ++z) {
        columnDensity += static_cast<double>(world.nH3p[z]) * static_cast<double>(world.dzcm[z]); // cm^-2
    }
    return columnDensity;
}

