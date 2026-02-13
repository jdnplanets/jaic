/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // H3pSpectrum.cpp

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstddef>
#include <iomanip>

#include "H3pSpectrum.hpp"


using namespace std;
namespace fs = filesystem;


// Global constants
const double k_B = 1.380649e-23;  // Boltzmann constant (J/K)
const double c = 2.99792458e8;   // Speed of light (m/s)
const double h = 6.62607015e-34; // Planck constant (J·s)

// H3pSpectrum class definition

// Constructor: Model initialization (read files etc.)
H3pSpectrum::H3pSpectrum() {

    // Here we perform initialization that is not specific to a given run

    cout << "Initializing H3pSpectrum..." << endl;
    cout << "Reading line data..." << endl;

    // Read data files
    read_line_data("h3p_line_list_neale_1996_subset.txt");

    num_points = static_cast<size_t>((lambda_max - lambda_min) / d_lambda) + 1;

    cout << "H3pSpectrum initialized." << endl;
}

// // generate() computes the H3+ spectrum for a given set of input parameters T and n.
void H3pSpectrum::generate(double T, double n, vector<double>& lambda, vector<double>& spec) {

    // cout << "" << endl;
    // cout << "****************************************************************" << endl;
    // cout << "*                   H3+ NIR Spectrum Model                     *" << endl;
    // cout << "****************************************************************" << endl;
    // cout << "Input parameters:" << endl;
    // cout << "Temperature: " << T << " K" << endl;
    // cout << "H3+ density: " << n << " m^-2" << endl;
    // cout << "****************************************************************" << endl;
    // cout << "" << endl;

    this->T = T;
    this->n = n;

    computeQ();

    emission_rates.clear();
    compute_line_emission_rates();

    
    // Now put everything on a grid and convolve with a Gaussian
    lambda_out.clear();
    spec_out.clear();

    // cout << "Generating spectrum on grid..." << endl;
    // grid_convolve_spectrum();
    grid_spectrum();

    // Write the spectrum to the output vectors
    for (size_t i = 0; i < num_points; ++i) {
        lambda[i] = lambda_out[i];
        spec[i]   = spec_out[i];
    }

    // for (size_t i = 0; i < num_points; i++) {
    //     cout << "  λ = " << lambda_out[i] << " μm, I = " << spec_out[i] << " W m^-2 sr^-1 μm^-1" << endl;
    // }

    // printf("Finished.\n");

}

// Function to write the final spectrum to a file
void H3pSpectrum::write_spectrum_to_file() {
    ofstream output(spec_outfile, ios::out | ios::trunc);
    output << "Wavelength (Å), Emission Rate" << endl;
    for (size_t i = 0; i < num_points; i++) {
        output << setw(8) << lambda_out[i] << " " << setw(12) << spec_out[i] << endl;
    }
    output.close();
}


// -----------------------------------------------------------
// Data reading functions
// -----------------------------------------------------------

// Function to read the line list data
void H3pSpectrum::read_line_data(const string& filename) {
    ifstream file(datadir+filename);
    if (!file) {
        cerr << "Error opening file: " << datadir+filename << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        if (line[0] == '#') continue;

        istringstream iss(line);
        int Ju;
        double wu, wl, A, g;

        if (!(iss >> Ju >> wu >> wl >> A >> g)) {
            cerr << "Error reading line in file: " << filename << endl;
            continue;
        }

        lines.push_back({Ju, wu, wl, A, g});
    }
    file.close();
}



// -----------------------------------------------------------
// Spectrum computation functions
// -----------------------------------------------------------



void H3pSpectrum::computeQ() {
    if (T >= 100 && T < 1800) {
        Q = -1.11391 + 0.0581076 * T + 0.000302967 * pow(T,2) - 2.83724e-7 * pow(T,3)
            + 2.31119e-10 * pow(T,4) - 7.15895e-14 * pow(T,5) + 1.00150e-17 * pow(T,6);
    } 
    else if (T >= 1800 && T < 5000) {
        Q = -22125.5 + 51.1539 * T - 0.0472256 * pow(T,2) + 2.26131e-5 * pow(T,3)
            - 5.85307e-9 * pow(T,4) + 7.90879e-13 * pow(T,5) - 4.28349e-17 * pow(T,6);
    } 
    else if (T >= 5000 && T < 10000) {
        Q = -654293.0 + 617.630 * T - 0.237058 * pow(T,2) + 4.74466e-5 * pow(T,3)
            - 5.20566e-9 * pow(T,4) + 3.05824e-13 * pow(T,5) - 7.45152e-18 * pow(T,6);
    } 
    else {
        cerr << "Error: Partition constants out of range of temperature " << T << " K" << endl;
        Q = 0.0;
    }
}


// Function to compute the emission rates for each bound transition
void H3pSpectrum::compute_line_emission_rates(){
    emission_rates.clear();

    // cout << "Computing line emission rates..." << endl;
    // Compute emissions for each transition
    for (const auto& ln : lines) {
        
        double expnt   = ( -100 * ln.wu *  h * c) / (k_B * T);
        double emrate  = ln.g * (2 * ln.Ju + 1) * 100 * h *  c * ln.wl * ln.A  / (4 * M_PI * Q) * exp(expnt);

        // cout << "Line Ju=" << ln.Ju << ", Wavelength =" << 1e4/ln.wl << " μm, Emission Rate=" << emrate<< " s^-1" << endl;
        
        emission_rates.push_back({emrate * n, ln.wu, 1e4/ln.wl}); // Store emission rate, energy (cm^-1), wavelength (μm)

    }
}



//--------------------------------------------------------------------------------
// Functions to generate a spectrum on a grid from the lines
//--------------------------------------------------------------------------------

// Function to add a line to the spectrum
void H3pSpectrum::add_line(vector<double>& spectrum, double lambda_line, double emission_rate, 
                double lambda_min, double lambda_max, double d_lambda) {
    // Find the closest wavelength index
    size_t index = static_cast<size_t>((lambda_line - lambda_min) / d_lambda);

    // Check that the index is within the grid bounds
    if (lambda_line >= lambda_min && lambda_line <= lambda_max) {
        spectrum[index] += emission_rate / d_lambda;  // Convert to spectral density
        // cout << "Added line at " << lambda_line << " μm to index " << index << " with emission rate " << emission_rate << endl;
    }
}

// Function to convolve the spectrum with a Gaussian kernel
vector<double> H3pSpectrum::convolve_with_gaussian(vector<double> spectrum, double d_lambda, double fwhm) {
    double sigma = fwhm / (2 * sqrt(2 * log(2)));  // Convert FWHM to sigma
    int kernel_half_width = static_cast<int>(3 * sigma / d_lambda);  // Kernel size ~ 3σ

    // Create and normalise the Gaussian kernel
    vector<double> gaussian_kernel(2 * kernel_half_width + 1);
    double kernel_sum = 0.0;
    for (int k = -kernel_half_width; k <= kernel_half_width; k++) {
        double weight = exp(-pow(k * d_lambda, 2) / (2 * pow(sigma, 2)));
        gaussian_kernel[k + kernel_half_width] = weight;
        kernel_sum += weight;
    }

    // normalise the kernel (correct area normalisation)
    for (double& weight : gaussian_kernel) {
        weight /= kernel_sum;
    }

    // Convolve the spectrum
    vector<double> convolved_spectrum(spectrum.size(), 0.0);
    for (size_t i = 0; i < spectrum.size(); i++) {
        for (int k = -kernel_half_width; k <= kernel_half_width; k++) {
            std::ptrdiff_t j = static_cast<std::ptrdiff_t>(i) + k;
            if (j >= 0 && j < static_cast<std::ptrdiff_t>(spectrum.size())) {
                convolved_spectrum[i] += spectrum[static_cast<std::size_t>(j)] * gaussian_kernel[k + kernel_half_width];
            }
        }
    }

    return convolved_spectrum;
}


void H3pSpectrum::grid_convolve_spectrum() {

    vector<double> spectrum(num_points, 0.0);
    vector<double> lambda(num_points, 0.0);
    // cout << "Creating wavelength grid from " << lambda_min << " μm to " << lambda_max << " μm with " << num_points << " points." << endl;
    for (size_t i = 0; i < num_points; i++) {
        lambda[i] = lambda_min + i * d_lambda;
    }

    // Add each emission line to the spectrum using a simple line or Doppler broadening
    for (const auto& line : emission_rates) {
        if (line.wavelength >= lambda_min && line.wavelength <= lambda_max) {
            // add_doppler_broadened_line(spectrum, lambda_line, line.emission_rate, T, lambda_min, lambda_max, d_lambda);
            add_line(spectrum, line.wavelength, line.emission_rate, lambda_min, lambda_max, d_lambda);
        }
    }

    // Convolve the spectrum with the Gaussian kernel of specified FWHM
    spectrum = convolve_with_gaussian(spectrum, d_lambda, fwhm);

    for (size_t i = 0; i < num_points; i++) {
        spec_out.push_back(spectrum[i]);
        lambda_out.push_back(lambda[i]);
    }
}

void H3pSpectrum::grid_spectrum() {

    vector<double> spectrum(num_points, 0.0);
    vector<double> lambda(num_points, 0.0);
    // cout << "Creating wavelength grid from " << lambda_min << " μm to " << lambda_max << " μm with " << num_points << " points." << endl;
    for (size_t i = 0; i < num_points; i++) {
        lambda[i] = lambda_min + i * d_lambda;
    }

    // Add each emission line to the spectrum using a simple line or Doppler broadening
    for (const auto& line : emission_rates) {
        if (line.wavelength >= lambda_min && line.wavelength <= lambda_max) {
            // add_doppler_broadened_line(spectrum, lambda_line, line.emission_rate, T, lambda_min, lambda_max, d_lambda);
            add_line(spectrum, line.wavelength, line.emission_rate, lambda_min, lambda_max, d_lambda);
        }
    }

    for (size_t i = 0; i < num_points; i++) {
        spec_out.push_back(spectrum[i]);
        lambda_out.push_back(lambda[i]);
    }
}



