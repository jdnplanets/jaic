/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file implements the H3pColumn class, which models the vertical column of Jupiter's ionosphere 
to compute the H3+ emission spectrum and related quantities. It uses the H3pSpectrum class to 
compute the spectrum based on input parameters such as temperature and H3+ density profiles.
The class also includes methods for applying non-LTE scaling factors, computing volume emission rates, 
and writing output to files. 
*/


 // H3pColumn.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <limits>

#include "H3pColumn.hpp"
#include "FileIO.hpp"


using namespace std;


// Constructor: Model initialization (read files, set up normalisation factors, etc.)
H3pColumn::H3pColumn(Params params, std::shared_ptr<Source> src_, const World1D &world_): 
                sp(params), src(src_), world(world_) {

    // Here we perform initialization that is not specific to a given run
    // For example, initialize g_J and E_J.
    cout << "Initializing H3pColumn..." << endl;

    Z = world.Z;    // Altitude grid in m
    nz = world.nz;  // Number of altitude bins
    dz = world.dz;  // Size of each bin in m
    T = world.T;    // Temperature profile in K
    n_H3p.assign(world.nH3p.begin(), world.nH3p.end()); // Neutral H3p density in cm^-3 (need to do this since nH3p is a vector<float> in World1D)
    // Convert n_H3p from cm^-3 to m^-3
    for (auto &val : n_H3p) val *= 1e6;

    // Initialise the output variables
    nw = H3pspec.num_points;
    specout.assign(nw, 0.0);
    ver.assign(nz, 0.0);
    q01ver.assign(nz, 0.0);

    // output file path
    outdir += sp.runid+"/";
    // Check if the directory exists, and create it if it doesn't
    if (!std::filesystem::exists(outdir)) {
        if (!std::filesystem::create_directories(outdir)) {
            throw std::runtime_error("Failed to create output directory: " + outdir);
        }
    }

    lambda.resize(nw, 0.0);
    H3pspecvz.clear();

    // compute non-LTE scaling factors
    nonLTEs.resize(nz, 1.0);
    for (size_t z = 0; z < nz; z++) {
        nonLTEs[z] = getNonLTEScalingFactor(z);
    }
}

double H3pColumn::getNonLTEScalingFactor(size_t z) {
    // Use the lookup table in H3pColumn to bilinearly interpolate a scaling
    // factor based on temperature and H2 density (world.nH2 is in cm^-3).
    if (z >= static_cast<size_t>(nz)) return 1.0;
    double Tval = T[z];
    double nH2_val = 0.0;
    if (static_cast<size_t>(world.nH2.size()) > z) {
        nH2_val = world.nH2[z] * 1e6; // m^-3
    } else {
        return 1.0;
    }

    const auto &temps = H3pColumn::nonLTE_temps;
    const auto &dens = H3pColumn::nonLTE_h2dens;
    const auto &scal = H3pColumn::nonLTE_scalings;

    // Clamp or locate temperature index
    size_t ti = 0;
    if (Tval <= temps.front()) {
        ti = 0;
    } else if (Tval >= temps.back()) {
        ti = temps.size() - 2;
        if (temps.size() < 2) ti = 0;
    } else {
        for (size_t i = 0; i + 1 < temps.size(); ++i) {
            if (Tval >= temps[i] && Tval <= temps[i+1]) { ti = i; break; }
        }
    }

    // Use log density for interpolation across many decades
    double logn = std::log10(std::max(nH2_val, 1e-30));
    std::vector<double> logdens(dens.size());
    for (size_t j = 0; j < dens.size(); ++j) logdens[j] = std::log10(dens[j]);

    size_t dj = 0;
    if (logn <= logdens.front()) {
        dj = 0;
    } else if (logn >= logdens.back()) {
        dj = logdens.size() - 2;
        if (logdens.size() < 2) dj = 0;
    } else {
        for (size_t j = 0; j + 1 < logdens.size(); ++j) {
            if (logn >= logdens[j] && logn <= logdens[j+1]) { dj = j; break; }
        }
    }

    // Fractions along each axis
    double t = 0.0;
    if (temps.size() >= 2) {
        double t0 = temps[ti];
        double t1 = temps[std::min(ti + 1, temps.size() - 1)];
        if (t1 > t0) t = (Tval - t0) / (t1 - t0);
    }

    double u = 0.0;
    if (logdens.size() >= 2) {
        double u0 = logdens[dj];
        double u1 = logdens[std::min(dj + 1, logdens.size() - 1)];
        if (u1 > u0) u = (logn - u0) / (u1 - u0);
    }

    // Fetch four neighbors
    double v00 = scal[ti][dj];
    double v10 = scal[std::min(ti + 1, scal.size() - 1)][dj];
    double v01 = scal[ti][std::min(dj + 1, scal[0].size() - 1)];
    double v11 = scal[std::min(ti + 1, scal.size() - 1)][std::min(dj + 1, scal[0].size() - 1)];

    // Bilinear interpolation
    double value = 
        (1 - t) * (1 - u) * v00 + t * (1 - u) * v10 + (1 - t) * u * v01 + t * u * v11;
    return value;
}

void H3pColumn::computeSpectraVsAltitude() {

    for (size_t z = 0; z < nz; z++) {

        // if (z != 40) continue;
        cout << "\rComputing spectrum for z = " << Z[z]/1e3 << " km" << flush;
        vector<double> spec;
        spec.assign(nw, 0.0); // Initialize the spectrum for this altitude
        H3pspec.generate(T[z], n_H3p[z], lambda, spec); // spec in W m^-2 sr^-1 μm^-1 for 1 m slab
        for (size_t i = 0; i < nw; i++) 
            spec[i] *= nonLTEs[z]; // Apply non-LTE scaling factor uniformly to all wavelengths - a big approximation!
        

        H3pspecvz.push_back(spec);
    }
    cout << endl; // New line after the loop
}


// Function to compute the volume emission rate in W m-3 sr-1 versus altitude
void H3pColumn::computeVERvsAltitude() {

        for (size_t z = 0; z < nz; z++) {
            double sum = 0.0;
            for (size_t w = 0; w < H3pspecvz[z].size(); w++) {
                sum += H3pspecvz[z][w] * H3pspec.d_lambda;
            }
            ver[z] = sum;
        }

}


// Function to compute the volume emission rate of the Q(0,1-) line in W m-3 sr-1 versus altitude
void H3pColumn::computeQ01VERvsAltitude() {

        for (size_t z = 0; z < nz; z++) {
            double sum = 0.0;
            for (size_t w = 0; w < H3pspecvz[z].size(); w++) {
                if (lambda[w] >= 3.948 && lambda[w] <= 3.958) // Q(0,1-) line
                    sum += H3pspecvz[z][w] * H3pspec.d_lambda;
            }
            q01ver[z] = sum;
        }

}

void H3pColumn::computeEmergentSpectrum(){

    for (size_t w = 0; w < nw; w++) {
        for (size_t z = 0; z < nz; z++) {
            // The contribution from this layer.
            specout[w] += H3pspecvz[z][w] * dz; // W m^-2 sr^-1 μm^-1
        }
    }
    specout = H3pspec.convolve_with_gaussian(specout, H3pspec.d_lambda, H3pspec.fwhm); // Convolve the spectrum with a Gaussian kernel
}




// Function to compute the total radiance in W m^-2 sr^-1
double H3pColumn::getTotalRadiance() {
    double total = 0.0;
    for (size_t w = 0; w < nw; w++) {
        total += specout[w] * H3pspec.d_lambda; // W m^-2 sr^-1
    }
    if (total == 0.0) {
        cerr << "Warning: Total radiance is zero, check the model parameters." << endl;
        return 0.0;
    }

    return total; 
}
// Function to compute the radiance in W m^-2 sr^-1 in the JWST NIRCam F335M filter band (3.18 to 3.54 μm)
double H3pColumn::getF335MRadiance() {
    double total = 0.0;
    for (size_t w = 0; w < nw; w++) {
        if (lambda[w] >= 3.18 && lambda[w] <= 3.54) {
            total += specout[w] * H3pspec.d_lambda; // W m^-2 sr^-1
        }
    }
    if (total == 0.0) {
        cerr << "Warning: F335M radiance is zero, check the model parameters." << endl;
        return 0.0;
    }

    return total; 
}


int H3pColumn::getPeakzIndex() {
    if (nz == 0) return -1;
    double max_ver = std::numeric_limits<double>::lowest();
    int peak_index = -1;
    for (size_t z = 0; z < nz; ++z) {
        if (ver[z] > max_ver) {
            max_ver = ver[z];
            peak_index = static_cast<int>(z);
        }
    }
    return peak_index;
}



void H3pColumn::writeVERtoFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "H3pver_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing VER to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H3+ VER"},
        {"layout", "rows=z_index (0..nz-1)"},
        {"nz", std::to_string(nz)},
        {"Zgrid_type", "linear"},
        {"Zmin (m)", std::to_string(Z.empty() ? 0.0 : Z.front())},
        {"Zmax (m)", std::to_string(Z.empty() ? 0.0 : Z.back())},
        {"Z units", "km"},
        {"VER units", " W m-3 sr-1"} // adjust if different
    };

    const std::vector<std::string> cols = {"Z [km]", "VER [W m-3 sr-1]"};

    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H3pColumn model output: H3+ volume emission rate (VER)",
        meta,
        cols,
        static_cast<int>(nz),
        [&](int i, std::ostream& os, int w) {
            os << std::right
               << std::setw(w) << static_cast<int>(Z[i] / 1e3) << " "
               << std::setw(w) << ver[i];
        },
        colw,
        precision
    );

    if (!ok) {
        std::cerr << "Failed to write H3p VER to " << filename << "\n";
    }
}


// void H3pColumn::writeQ01VERtoFile() {
//     stringstream ss;
//     if (!src) { cerr << "Error: src is null" << endl;return; }
//     string suf = src->label(sp);
//     ss << outdir << "H3pq01ver_" << suf << ".dat";
//     string filename = ss.str();
//     ofstream outFile(filename, ios::trunc);
//     if (!outFile) {
//         cerr << "Error opening file for writing: " + filename << endl;
//         return;
//     }
//     cout << "Writing VER to file: " << filename << endl;

//     for (size_t i = 0; i < nz; ++i) {
//         outFile << Z[i]/1e3 << " " << q01ver[i] << endl;
//     }


//     outFile.close();
// }


void H3pColumn::writeQ01VERtoFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "H3pq01ver_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing Q(0,1-) VER to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H3+ Q(0,1-) VER"},
        {"layout", "rows=z_index (0..nz-1)"},
        {"nz", std::to_string(nz)},
        {"Zgrid_type", "linear"},
        {"Zmin (m)", std::to_string(Z.empty() ? 0.0 : Z.front())},
        {"Zmax (m)", std::to_string(Z.empty() ? 0.0 : Z.back())},
        {"Z units", "km"},
        {"VER units", " W m-3 sr-1"} // adjust if different
    };

    const std::vector<std::string> cols = {"Z [km]", "Q(0,1-) VER [W m-3 sr-1]"};

    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H3pColumn model output: H3+ Q(0,1-) line volume emission rate (VER)",
        meta,
        cols,
        static_cast<int>(nz),
        [&](int i, std::ostream& os, int w) {
            os << std::right
               << std::setw(w) << static_cast<int>(Z[i] / 1e3) << " "
               << std::setw(w) << q01ver[i];
        },
        colw,
        precision
    );

    if (!ok) {
        std::cerr << "Failed to write H3p Q(0,1-) VER to " << filename << "\n";
    }
}



// void H3pColumn::writeEmergentSpectrumToFile() {
//     stringstream ss;
//     if (!src) { cerr << "Error: src is null" << endl;return; }
//     string suf = src->label(sp);
//     ss << outdir << "H3pspecout_" << suf << ".dat";
//     string filename = ss.str();
//     ofstream outFile(filename, ios::trunc);
//     if (!outFile) {
//         cerr << "Error opening file for writing: " + filename << endl;
//         return;
//     }
//     cout << "Writing emergent spectrum to file: " << filename << endl;

//     for (size_t i = 0; i < nw; ++i) {
//         outFile << lambda[i] << " " << specout[i] << endl;
//     }

//     outFile.close();
// }


void H3pColumn::writeEmergentSpectrumToFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "H3pspecout_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing emergent spectrum to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H3+ emergent spectrum"},
        {"layout", "rows=w_index (0..nw-1)"},
        {"nw", std::to_string(nw)},
        {"wavelength units", "μm"},
        {"radiance units", "W m^-2 sr^-1 μm^-1"} 
    };

    const std::vector<std::string> cols = {"wavelength [μm]", "spectral radiance [W m^-2 sr^-1 μm^-1]"};

    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H3pColumn model output: H3+ emergent spectrum",
        meta,
        cols,
        static_cast<int>(nw),
        [&](int i, std::ostream& os, int w) {
            os << std::right
               << std::setw(w) << lambda[i]  << " "
               << std::setw(w) << specout[i];
        },
        colw,
        precision
    );

    if (!ok) {
        std::cerr << "Failed to write H3p emergent spectrum to " << filename << "\n";
    }
}
