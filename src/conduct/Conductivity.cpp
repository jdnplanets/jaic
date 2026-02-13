/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // Conductivity.cpp

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "World1D.hpp"
#include "Conductivity.hpp"



#include "constants.h"
using namespace std;
using namespace constants;

Conductivity::Conductivity(Params params_, 
                           std::shared_ptr<Source> src_, 
                           World1D &world_)
    : sp(params_), src(src_), world(world_){
    // Constructor initializes the world model, ionization data, and ion

    sigmaP_H3p.resize(world.nz, 0.0); // Initialize conductivity vector
    sigmaP_CH5p.resize(world.nz, 0.0);
    sigmaP_C3Hnp.resize(world.nz, 0.0);
    sigmaP_e.resize(world.nz, 0.0);
    sigmaP.resize(world.nz, 0.0); // Initialize conductivity vector

        // output file path
    outdir += sp.runid+"/";
    // Check if the directory exists, and create it if it doesn't
    if (!std::filesystem::exists(outdir)) {
        if (!std::filesystem::create_directories(outdir)) {
            throw std::runtime_error("Failed to create output directory: " + outdir);
        }
    }
}


void Conductivity::computeConductivity() {

    // Compute conductivity based on the ion densities and species properties.

    for (size_t z = 0; z < world.nz; ++z) {

        double T = world.T[z]; // Temperature in K
        double B = world.B; // Magnetic field strength in Tesla
        double n_e = world.ne[z] * 1e6; 
        double n_H2 = world.nH2[z] * 1e6;      


        // double sigma_en = pi * (H2.r + Electron.r) * (H2.r + Electron.r);       // electron-neutral collision cross-section in m^2
        double m_en = m_e * H2.m / (m_e + H2.m);                                // electron-neutral reduced mass in kg
        // double T_en = (m_e * T + H2.m * T) / (m_e + H2.m);                      // average temperature in K
        double Omega_ce = qe * B / (m_e);                                       // electron cyclotron frequency in rad/s
        double nu_en = 2.605e-10 * world.nH2[z] * pow(T, 0.6159);               // Given in Koskinen et al. 2010 n_H2 in cm^-3
        sigmaP_e[z] = n_e * qe * qe * 
                    (1 / (m_e* nu_en) * (nu_en * nu_en /(nu_en * nu_en + Omega_ce * Omega_ce)));       
        // cout << scientific << sigmaP_e[z] << " " << world.nH2[z] << endl;

        double n_i = world.nH3p[z] * 1e6;                                       // Ion density in m^-3
        double m_in = H3p.m * H2.m / (H3p.m + H2.m);                            // ion-neutral reduced mass in kg
        double Omega_ci = H3p.q * B / (H3p.m);                                  // ion cyclotron frequency in rad/s     
        double nu_in = 2.6e-15 * sqrt(0.82 / (m_in/m_p)) * n_H2;
        sigmaP_H3p[z] = n_i * qe * qe * 
                    (1 / (H3p.m * nu_in) * (nu_in * nu_in /(nu_in * nu_in + Omega_ci * Omega_ci)));      


        n_i = world.nCH5p[z] * 1e6;                                       // Ion density in m^-3
        m_in = CH5p.m * H2.m / (CH5p.m + H2.m);                            // ion-neutral reduced mass in kg
        Omega_ci = CH5p.q * B / (CH5p.m);                                  // ion cyclotron frequency in rad/s
        nu_in = 2.6e-15 * sqrt(0.82 / (m_in/m_p)) * n_H2;
        sigmaP_CH5p[z] = n_i * qe * qe * 
                    (1 / (CH5p.m * nu_in) * (nu_in * nu_in /(nu_in * nu_in + Omega_ci * Omega_ci)));
 
        n_i = world.nC3Hnp[z] * 1e6;                             // Ion density in m^-3
        m_in = C3Hnp.m * H2.m / (C3Hnp.m + H2.m);                            // ion-neutral reduced mass in kg
        Omega_ci = C3Hnp.q * B / (C3Hnp.m);                                  // ion cyclotron frequency in rad/s
        nu_in = 2.6e-15 * sqrt(0.82 / (m_in/m_p)) * n_H2;
        sigmaP_C3Hnp[z] = n_i * qe * qe * 
                    (1 / (C3Hnp.m * nu_in) * (nu_in * nu_in /(nu_in * nu_in + Omega_ci * Omega_ci)));
 
        // Total Pedersen conductivity
        sigmaP[z] = sigmaP_e[z] + sigmaP_H3p[z] + sigmaP_CH5p[z] + sigmaP_C3Hnp[z];
    }
}


void Conductivity::writeConductivityToFile() const {
    stringstream ss;
    if (!src) { cerr << "Error: src is null" << endl;return; }
    string suf = src->label(sp);
    ss << outdir << "sigmap_" << suf << ".dat";
    string filename = ss.str();
    ofstream outFile(filename, ios::trunc);
    if (!outFile) {
        cerr << "Error opening file for writing: " + filename << endl;
        return;
    }

    outFile << "# Altitude (m), H3+ Conductivity (S/m), CH5+ Conductivity (S/m), C3Hn+ Conductivity (S/m), Electron Conductivity (S/m), Total Conductivity (S/m)\n";
    for (size_t z = 0; z < world.nz; ++z) {
        outFile << world.Z[z] << "\t "
                << sigmaP_H3p[z] << "\t "
                << sigmaP_CH5p[z] << "\t "
                << sigmaP_C3Hnp[z] << "\t "
                << sigmaP_e[z] << "\t "
                << sigmaP[z] << "\n";
    }
    outFile.close();
    cout << "Conductivity data written to " << filename << endl;
}


double Conductivity::getConductance() {
    // Integrate conductivity over altitude to get conductance (S)
    double SigmaP = 0.0;
    for (size_t z = 0; z < world.nz - 1; ++z) {
        double dz = world.Z[z+1] - world.Z[z]; // in m
        SigmaP += 0.5 * (sigmaP[z] + sigmaP[z+1]) * dz; // Trapezoidal rule
    }
    return SigmaP;
}