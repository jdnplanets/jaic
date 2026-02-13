/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
This file inmplements the H2Column class, which computes the H2 emission spectrum from a column of Jupiter's 
atmosphere given excitation rates and atmospheric composition. It includes methods to read excitation rates, 
compute volume emission rates, and calculate the emergent spectrum accounting for absorption by hydrocarbons. 
The class relies on the World1D model for atmospheric properties and a Source for excitation rates. 
*/

 // H2Column.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <filesystem>

#include "H2Column.hpp"
#include "Egrid.hpp"
#include "FileIO.hpp"


using namespace std;


// Constructor: Model initialization (read files, set up normalisation factors, etc.)
H2Column::H2Column(Params params, std::shared_ptr<Source> src_, const World1D &world_, const bool hydrocarbon_absorption): 
                sp(params), src(src_), world(world_), absorb(hydrocarbon_absorption) {

    // Here we perform initialization that is not specific to a given run
    // For example, initialize g_J and E_J.

    cout << "Initializing H2Column..." << endl;


    Z = world.Z;                // Altitude grid in m
    nz = world.nz;              // Number of altitude bins
    dz = world.dz;              // Size of each bin in m
    dzcm = world.dzcm;          // Size of each bin in cm
    T = world.T;                // Temperature profile in K
    n_H2.assign(world.nH2.begin(), world.nH2.end()); // Neutral H2 density in cm^-3 (need to do this since nH2 is a vector<float> in World1D)
    n_CH4 = world.nCH4;         // Methane density in cm^-3
    n_C2H2 = world.nC2H2;       // Acetylene density in cm^-3
    n_C2H4 = world.nC2H4;       // Ethylene density in cm^-3
    n_C2H6 = world.nC2H6;       // Ethane density in cm^-3

    col_CH4 = vector<double>(nz, 0.0);  // Column density of CH4 from each altitude to the top in cm^-2
    col_C2H2 = vector<double>(nz, 0.0); // Column density of C2H2 from each altitude to the top in cm^-2
    col_C2H4 = vector<double>(nz, 0.0); // Column density of C2H4 from each altitude to the top in cm^-2
    col_C2H6 = vector<double>(nz, 0.0); // Column density of C2H6 from each altitude to the top in cm^-2

    // Initialise the output variables
    nw = h2spec.num_points;

    // Resize vectors to the correct sizes (reserve alone is insufficient)
    lambda.assign(nw, 0.0);
    specout.assign(nw, 0.0);
    h2specvz.assign(nz, vector<double>(nw, 0.0));
    ver.assign(nz, 0.0);
    ver_eV.assign(nz, 0.0);

    // output file path
    outdir += sp.runid+"/";
    // Check if the directory exists, and create it if it doesn't
    if (!std::filesystem::exists(outdir)) {
        if (!std::filesystem::create_directories(outdir)) {
            throw std::runtime_error("Failed to create output directory: " + outdir);
        }
    }


    if (absorb) {

    // Precompute column density profiles (integrated from each altitude to the top)
    // We compute these using the trapezoidal rule, working downward (from the top).
        // At the top (last grid point), assume zero column density above.
        // // Integrate from bottom (lower altitudes) upward.
        CnHn_xs_data = readCnHnXSectData("interpolated_hydrocarbon_photon_xsects.csv");

        for (int i = nz - 2; i >= 0; --i) {
            col_CH4[i]  = col_CH4[i+1]  + 0.5 * (n_CH4[i]  + n_CH4[i+1])  * dzcm;
            col_C2H2[i] = col_C2H2[i+1] + 0.5 * (n_C2H2[i] + n_C2H2[i+1]) * dzcm;
            col_C2H4[i] = col_C2H4[i+1] + 0.5 * (n_C2H4[i] + n_C2H4[i+1]) * dzcm;
            col_C2H6[i] = col_C2H6[i+1] + 0.5 * (n_C2H6[i] + n_C2H6[i+1]) * dzcm;
        }
    }

    // Build default energy grid from simparams defaults
    {
        SimParams spSim; // use parameters from simparams.h
        energyGrid = Egrid::fromEminEmax(spSim.nbinsE, static_cast<double>(spSim.E0), static_cast<double>(spSim.E0) * std::pow(10.0, static_cast<double>(spSim.decades)));
    }
}

 void H2Column::setExcitationRates(vector<vector<double>>& RB, vector<vector<double>>& RC, vector<vector<double>>& RE, const double F) {
        R_B = RB;
        R_C = RC;
        R_E = RE;

        // If the excitation-energy axis doesn't match our current energy grid, adjust it now.
        if (!R_B.empty()) {
            size_t nbinsE_from_rates = R_B[0].size();
            if (static_cast<int>(nbinsE_from_rates) != energyGrid.nbins) {
                // Rebuild energyGrid to match the rates (use SimParams defaults for Emin/decades)
                SimParams spSim;
                energyGrid = Egrid::fromEminEmax(static_cast<int>(nbinsE_from_rates), static_cast<double>(spSim.E0), static_cast<double>(spSim.E0) * std::pow(10.0, static_cast<double>(spSim.decades)));
            }
        }

        if (F != 1.0) {
            for (size_t i = 0; i < R_B.size(); ++i) {
                for (size_t j = 0; j < R_B[i].size(); ++j) {
                    R_B[i][j] *= F;
                    R_C[i][j] *= F;
                    R_E[i][j] *= F;
                }
            }
        }
        
    }

// Function to read excitation rates from a file
void H2Column::readExcitationRates(const Params &params,
                                  const std::shared_ptr<Source> &src,
                                  double F)
{
    const std::string basepath = "out/precip/" + params.runid + "/";
    const std::string suf = src->label(sp);
    const std::string filename = basepath + "exrates_" + suf + ".dat";

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open excitation rates file: " + filename);
    }

    const int nbinsE = energyGrid.nbins;
    const int nz_expected = nz;

    R_B.assign(nz_expected, std::vector<double>(nbinsE, 0.0));
    R_C.assign(nz_expected, std::vector<double>(nbinsE, 0.0));
    R_E.assign(nz_expected, std::vector<double>(nbinsE, 0.0));

    std::string line;
    int z = 0;

    while (std::getline(infile, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;

        if (z >= nz_expected) {
            // Extra data lines
            break;
        }

        std::istringstream iss(line);
        for (int e = 0; e < nbinsE; ++e) {
            double b, c, eF;
            if (!(iss >> b >> c >> eF)) {
                throw std::runtime_error(
                    "Malformed data row (z=" + std::to_string(z) + ") in " + filename +
                    " (expected " + std::to_string(3 * nbinsE) + " numbers)."
                );
            }
            R_B[z][e] = b;
            R_C[z][e] = c;
            R_E[z][e] = eF;

        }

        ++z;
    }

    if (z != nz_expected) {
        throw std::runtime_error(
            "Excitation rates file has " + std::to_string(z) +
            " data rows but expected nz=" + std::to_string(nz_expected) +
            " in " + filename
        );
    }

    // Scale by number flux
    for (int iz = 0; iz < nz_expected; ++iz) {
        for (int ie = 0; ie < nbinsE; ++ie) {
            R_B[iz][ie] *= F;
            R_C[iz][ie] *= F;
            R_E[iz][ie] *= F;
        }
    }

    std::cout << "Excitation rates read successfully from " << filename << "\n";
}


// Function to compute the spectrum vs altitude
// Spectrum units are photons cm^-2 s^-1 radiating into 4π sr
void H2Column::computeSpectraVsAltitude() {

    for (size_t z = 0; z < nz; z++) {

        cout << "\rComputing spectrum for z = " << Z[z]/1e3 << " km" << flush;
        vector<double> spec;
        spec.assign(nw, 0.0); // Initialize the spectrum for this altitude
        size_t nE = R_B[z].size(); // Number of energy levels
        // Ensure the rates' energy axis matches our energyGrid; warn if not.
        if (static_cast<int>(nE) != energyGrid.nbins) {
            cerr << "Warning: excitation rates energy axis (" << nE << ") does not match energyGrid.nbins (" << energyGrid.nbins << ").\n";
        }
        for (size_t e = 0; e < nE; e++) {
            if (R_B[z][e] < 1e-2) {
                continue; // Skip if B excitation is negligible
            }
            vector<double> specE;
            specE.assign(nw, 0.0);
            double E = energyGrid.eToE(static_cast<int>(e));
            h2spec.generate(T[z], E, R_B[z][e], R_C[z][e], R_E[z][e], lambda, specE);

            // Add the contribution from this energy level to the spectrum
           for (size_t i = 0; i < nw; i++) {
                spec[i] += specE[i];
            }
        }
        h2specvz[z] = spec;
    }
    cout << endl; // New line after the loop
}




// Function to read hydrocarbon photon cross section data from CSV file
// Cross sections are in units of cm^2/molecule
vector<H2Column::CnHnXSectData> H2Column::readCnHnXSectData(const string& filename) {
    vector<CnHnXSectData> data;
    ifstream file("data/h2spec/"+filename);
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    // Skip the header line
    getline(file, line);

    while (getline(file, line)) {
        istringstream iss(line);
        string value;
        CnHnXSectData point;
        
        getline(iss, value, ',');
        point.lambda = stod(value);

        getline(iss, value, ',');
        point.CH4 = stod(value);

        getline(iss, value, ',');
        point.C2H2 = stod(value);

        getline(iss, value, ',');
        point.C2H4 = stod(value);

        getline(iss, value, ',');
        point.C2H6 = stod(value);

        data.push_back(point);
    }

    return data;
}



// Given a wavelength (in Å) and a vector of cross‐section data (assumed sorted
// in ascending order by lambda), return the cross sections (for all species) 
// from the grid point nearest to the requested wavelength.
H2Column::HydrocarbonCrossSections 
            H2Column::getCrossSections(double wavelength,
                                            const vector<CnHnXSectData>& xsData)
{
    auto it = lower_bound(xsData.begin(), xsData.end(), wavelength,
        [](const CnHnXSectData& data, double w) { return data.lambda < w; });
    if(it == xsData.end()){
        const CnHnXSectData& d = xsData.back();
        return { d.CH4, d.C2H2, d.C2H4, d.C2H6 };
    }
    if(it == xsData.begin()){
        const CnHnXSectData& d = *it;
        return { d.CH4, d.C2H2, d.C2H4, d.C2H6 };
    }
    auto prev = it - 1;
    if( abs(it->lambda - wavelength) < abs(wavelength - prev->lambda) ){
        const CnHnXSectData& d = *it;
        return { d.CH4, d.C2H2, d.C2H4, d.C2H6 };
    } else {
        const CnHnXSectData& d = *prev;
        return { d.CH4, d.C2H2, d.C2H4, d.C2H6 };
    }
}


// Function to compute the volume emission rate in photons cm^-3 s^-1 radiating into 4π sr versus altitude
void H2Column::computeVERvsAltitude() {

        for (size_t z = 0; z < nz; z++) {
            double sum = 0.0;
            for (size_t w = 0; w < h2specvz[z].size(); w++) {
                sum += h2specvz[z][w];
            }
            ver[z] = sum;
        }

}

// Function to compute the volume emission rate in eV cm^-3 s^-1 versus altitude
void H2Column::computeVEReVvsAltitude() {

        for (size_t z = 0; z < nz; z++) {
            double sum = 0.0;
            for (size_t w = 0; w < h2specvz[z].size(); w++) {
                sum += h2specvz[z][w] * constants::h * constants::c / (lambda[w] * 1e-10) / constants::eV; // Convert wavelength from Å to m
            }
            ver_eV[z] = sum;
        }
}


// Function to compute the emergent spectrum at the top of the atmosphere
// Units are photons cm^-2 s^-1 radiating into 4π sr
void H2Column::computeEmergentSpectrum(){

    for (size_t i = 0; i < nw; i++) {
        double w = lambda[i];
        // Look up cross sections (for all species) at this wavelength. cm^2/molecule
        HydrocarbonCrossSections xs;
        if (absorb) xs = getCrossSections(w, CnHn_xs_data);

        for (size_t z = 0; z < nz; z++) {
            double tau = 0.0;
            if (absorb){
                tau = ( xs.CH4  * col_CH4[z] +
                        xs.C2H2 * col_C2H2[z] +
                        xs.C2H4 * col_C2H4[z] +
                        xs.C2H6 * col_C2H6[z] )
                        / cos(obsang * M_PI / 180);
            }
            // The contribution from this layer.
            specout[i] += h2specvz[z][i] * exp(-tau) * dzcm; // Convert m to cm. ds?
        }
    }
    specout = h2spec.convolve_with_gaussian(specout, h2spec.d_lambda, h2spec.fwhm); // Convolve the spectrum with a Gaussian kernel
}



// Function to compute the total unabsorbed intensity in kRayleighs
double H2Column::getTotalUnabsorbedIntensity() {
    double total = 0.0;
    for (size_t z = 0; z < nz; z++) {
        total += ver[z] * dzcm;
    }
    if (total == 0.0) {
        cerr << "Warning: Total intensity is zero, check the model parameters." << endl;
        return 0.0;
    }
    totkr0 = total / 1e9; // Convert to kRayleighs (1 kR = 1e9 photons cm^-2 s^-1 emitting into 4π sr)
    return totkr0; 
}

// Function to compute the total observed intensity in kRayleighs
double H2Column::getTotalObservedIntensity() {
    double total = 0.0;
    for (size_t w = 0; w < nw; w++) {
        total += specout[w];
    }
    if (total == 0.0) {
        cerr << "Warning: Total intensity is zero, check the model parameters." << endl;
        return 0.0;
    }
    totkr = total / (1e9*cos(obsang * M_PI / 180)); // Convert to kRayleighs (1 kR = 1e9 photons cm^-2 s^-1)
    return totkr; 
}

double H2Column::getColourRatio() {
    // Calculate the total emission between 1230 and 1300 Å
    double cr2 = 0.0;
    for (size_t w = 0; w < nw; w++) {
        if (lambda[w] >= 1230.0 && lambda[w] <= 1300.0) {
            cr2 += specout[w];
        }
    }
    if (cr2 == 0.0) {
        cerr << "Warning: Total emission in the 1230-1300 Å range is zero, colour ratio cannot be computed." << endl;
        return 0.0;
    }
    // Calculate the total emission between 1550 and 1620 Å
    double cr1 = 0.0;
    for (size_t w = 0; w < nw; w++) {
        if (lambda[w] >= 1550.0 && lambda[w] <= 1620.0) {
            cr1 += specout[w];
        }
    }
    cr = cr1 / cr2; // Colour ratio

    return cr;
}




void H2Column::writeVERtoFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "h2ver_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing VER to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H2 VER"},
        {"layout", "rows=z_index (0..nz-1)"},
        {"nz", std::to_string(nz)},
        {"Zgrid_type", "linear"},
        {"Zmin (m)", std::to_string(Z.empty() ? 0.0 : Z.front())},
        {"Zmax (m)", std::to_string(Z.empty() ? 0.0 : Z.back())},
        {"Z_units", "km"},
        {"VER_units", "ph cm-3 s-1 radiating into 4π sr"}
    };

    const std::vector<std::string> cols = {"Z [km]", "VER [ph cm-3 s-1]"};

    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H2Column model output: volume emission rate (VER)",
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
        std::cerr << "Failed to write VER to " << filename << "\n";
    }
}



void H2Column::writeVEReVtoFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "h2vereV_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing VER (eV) to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H2 VER (eV)"},
        {"layout", "rows=z_index (0..nz-1)"},
        {"nz", std::to_string(nz)},
        {"Zgrid_type", "linear"},
        {"Zmin (m)", std::to_string(Z.empty() ? 0.0 : Z.front())},
        {"Zmax (m)", std::to_string(Z.empty() ? 0.0 : Z.back())},
        {"Z_units", "km"},
        {"VER_units", "eV cm-3 s-1 radiating into 4π sr"}
    };

    const std::vector<std::string> cols = {"Z [km]", "VER [eV cm-3 s-1]"};
    
    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H2Column model output: volume emission rate (VER)",
        meta,
        cols,
        static_cast<int>(nz),
        [&](int i, std::ostream& os, int w) {
            os << std::right
               << std::setw(w) << static_cast<int>(Z[i] / 1e3) << " "
               << std::setw(w) << ver_eV[i];
        },
        colw,
        precision
    );

    if (!ok) {
        std::cerr << "Failed to write VER (eV) to " << filename << "\n";
    }
}




void H2Column::writeEmergentSpectrumToFile() {
    if (!src) {
        std::cerr << "Error: src is null\n";
        return;
    }

    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "h2specout_" << suf << ".dat";
    const std::string filename = ss.str();

    std::cout << "Writing emergent spectrum to file: " << filename << "\n";

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "H2 emergent_spectrum"},
        {"layout", "rows=w_index (0..nw-1)"},
        {"nw", std::to_string(nw)},
        {"lambda_units", "Angstrom"},
        {"specout_units", "ph cm^-2 s^-1 radiating into 4π sr"}
    };

    const std::vector<std::string> cols = {"wavelength [A]", "flux [ph cm^-2 s^-1]"};

    const int colw = 14;
    const int precision = 6;

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "H2Column model output: H2 emergent spectrum",
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
        std::cerr << "Failed to write emergent spectrum to " << filename << "\n";
    }
}
