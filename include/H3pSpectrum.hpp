/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // H3pSpectrum.hpp

#ifndef H3pSPECTRUM_H
#define H3pSPECTRUM_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>


class H3pSpectrum {
public:
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    std::string outdir = jaicroot + "out/h3pspec/";
    std::string spec_outfile;  // Output file for the spectrum
    static inline const std::string  datadir = jaicroot + "data/h3pspec/";  // Directory containing the data files
    // std::string lines_outfile = "synth_H3p_lines.dat";  // Output file for the spectrum
    // std::string spec_outfile = "synth_H3p_spec.dat";  // Output file for the spectrum
    double threshold = 1e-8;     // Threshold for negligible emission rates
    double lambda_min = 2.0;     // Minimum wavelength in μm
    double lambda_max = 5.0;     // Maximum wavelength in μm
    double d_lambda = 0.0001;      // Wavelength resolution in μm
    double fwhm = 0.0045;           // FWHM of the Gaussian in μm
    size_t num_points; // Number of points in the output spectrum
    std::vector<double> spec_out; // Spectrum grid


    // Constructor
    H3pSpectrum();

    // Member functions

    // generate() performs the spectrum calculations using input parameters.
    // T: temperature; n: H3+ density
    // The resulting convolved spectrum is stored in convolved_spectrum.
    void generate(double T, double n, std::vector<double>& lambda, std::vector<double>& spec);
    // void write_lines_to_file();
    void write_spectrum_to_file();
    std::vector<double> convolve_with_gaussian(std::vector<double> spectrum, double d_lambda, double fwhm);

protected:
    
    double T; // Temperature
    double n; // H3+ density
    double Q; // Partition function

    // Structure for storing line data
    struct Line {
        int Ju;
        double wu, wl, A, g;
    };

    // Structure to hold emission data
    struct EmissionData {
        double emission_rate;
        double energy;  // Energy of the upper state in cm^-1
        double wavelength;  // Wavelength of the emission line in Å
    };

    std::vector<Line> lines; // Transitions
    std::vector<EmissionData> emission_rates; // Emission rates
    std::vector<double> lambda_out; // Wavelength grid

    void read_line_data(const std::string& filename);
    void computeQ();
    void compute_line_emission_rates();

    void add_line(std::vector<double>& spectrum, double lambda_line, double emission_rate, 
                  double lambda_min, double lambda_max, double d_lambda);
    // void add_doppler_broadened_line(std::vector<double>& spectrum, double lambda_line, double emission_rate, 
                                    // double T, double lambda_min, double lambda_max, double d_lambda); 
    
    void grid_convolve_spectrum();
    void grid_spectrum();


};
#endif // H3pSPECTRUM_H