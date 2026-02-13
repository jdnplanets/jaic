/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
 This file defines the H2Spectrum class, which calculates the H2 emission spectrum based on excitation rates
and transition probabilities. It reads in transition data, computes excitation rates, and generates the emission spectrum,
 which can be convolved with a Gaussian profile to simulate instrumental broadening. 
 */

 // H2Spectrum.hpp

#ifndef H2SPECTRUM_H
#define H2SPECTRUM_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>

#include "read_continuum_data.h"


class H2Spectrum {
public:

    // File paths
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    std::string outdir = jaicroot + "out/h2spec/";
    static inline const std::string  datadir = jaicroot + "data/h2spec/";


    int max_J_ef = 32;              // Maximum J level for the EF state
    double Athresh;                 // Threshold for negligible transition probabilities
    double threshold = 1e-12;       // Threshold for negligible emission rates
    bool redistribute_thick = true; // Flag to redistribute optically thick emission rates
    double lambda_min = 875.0;      // Minimum wavelength in Å (87.5 nm)
    double lambda_max = 1800.0;     // Maximum wavelength in Å (180 nm)
    double d_lambda = 0.01;         // Wavelength resolution in Å
    double fwhm = 0.13;             // FWHM of the Gaussian in Å
    std::size_t num_points;         // Number of points in the output spectrum
    std::vector<double> spec_out;   // Spectrum grid


    // Constructor
    H2Spectrum();

    // Member functions

    // generate() performs the spectrum calculations using input parameters.
    // T: temperature; R_B, R_C, and R_E: excitation rates for different states.
    void generate(double T, double E, double R_B, double R_C, double R_E,  std::vector<double>& lambda, std::vector<double>& spec);
    void write_lines_to_file();
    void write_spectrum_to_file();
    void print_colour_ratio();
    void grid_convolve_spectrum();
    void grid_spectrum();
    std::vector<double> convolve_with_gaussian(std::vector<double> spectrum, double d_lambda, double fwhm);


protected:

    double total_excitation_rate;
    double totH2;
    
    // Structure for storing transition data
    struct Transition {
        int vu, Ju, vl, Jl;
        double Abound, energy, termu, terml;
        char branch;
    };

    // Structure for storing Franck–Condon factors from Abgrall and Roueff
    struct FranckCondon {
        int vu, Ju, vl, Jl;
        double fcf, energy;
        char branch;
    };

    // Structure to hold emission data
    struct EmissionData {
        double emission_rate;
        double energy;  // Energy of the upper state in cm^-1
        double wavelength;  // Wavelength of the emission line in Å
    };

    // Maps for transitions and Franck–Condon factors
    std::map<std::string, std::vector<Transition>> transitions;
    std::map<std::string, std::vector<FranckCondon>> franck_condon_transitions; // Franck–Condon factors from Abgrall and Roueff
    std::map<std::string, std::vector<std::vector<double>>> franck_condon_factors; // Franck–Condon factors from Fantz and Wunderlich


    // Map to store total excitation rate per level (v', J')
    std::map<std::string, std::map<std::pair<int, int>, double>> normalisation_factors_ef; // Normalisation factors for the EF state
    std::map<std::string, std::map<std::pair<int, int>, double>> level_excitation_rates_ef; // Excitation rates for the EF state
    std::map<std::string, std::map<std::tuple<int, int, int, int>, double>> EF_B_transition_energies; // EF to B transition energies
    std::map<std::string, std::map<std::pair<int, int>, double>> normalisation_factors_ef_cascade; // Normalisation factors for the EF->B cascade
    std::map<std::string, std::map<std::pair<int, int>, double>> level_excitation_rates; // Excitation rates for the B and C states
    std::map<std::string, std::map<std::pair<int, int>, double>> normalisation_factors_bc_upper; // Normalisation factors for the upper levels of the B and C states
    std::map<std::string, std::map<std::pair<int, int>, double>> normalisation_factors_bc_lower; // Normalisation factors for the lower levels of the B and C states
    std::map<std::string, std::map<std::pair<int, int>, double>> level_Acont_sums; // Sum of A coefficients for continuum transitions from each level
    std::map<std::string, std::map<std::pair<int, int>, double>> level_Abound_sums; // Sum of A coefficients for bound transitions from each level
    std::map<std::tuple<std::string, int, int, int, int>, EmissionData> emission_rates; // Emission rates for each transition, keyed by (state, vu, Ju, vl, Jl)

    // Map to store the continuum data
    ContinuumMap cmap;
    std::map<std::string, std::map<std::pair<int, int>, double>> continuum;
    std::vector<double> continuum_spec;

    std::vector<double> g_J;        // Statistical weights for ground state
    std::vector<double> E_J;        // Corresponding energy levels (in cm^-1)
    int max_Jl;                     // Maximum J level for the ground state
    std::vector<int> J_levels;      // J levels for the ground state
    int max_J_B;                    // Maximum J level for the B state
    std::vector<double> lambda_out; // Wavelength grid

// Member function declarations
    void read_transition_data(const std::string& filename, std::string state);
    void filter_negligible_transitions(double threshold);
    void read_franck_condon_factors_fw(const std::string& filename, std::string state);
    void read_franck_condon_factors_ar(const std::string& filename, std::string state);
    void read_continuum_data(const std::string& filename);
    double honl_london(int J, char branch, std::string state);
    double Abgrall_SJJ(int Jl, int Ju);
    double Ajello_SJJ(int Jl, int Ju);
    double Liu_SJJ(int Jl, int Ju);
    double Liu_correction_factor(const Transition& tr, std::string state);
    double compute_partition_function(const std::vector<double>& g_J, const std::vector<double>& E_J, double temp);
    void sum_Aconts();
    void sum_Abounds();
    void compute_EF_B_transition_energies();
    void excitation_XEF(double E, double R_E=0, std::vector<double> P_J={}, bool norm=false);
    void compute_normalisation_factors_excitation_XEF(double E);
    void cascade_EFB(bool norm=false);
    void compute_normalisation_factors_cascade_EFB();
    void excitation_XBC(double E, double R_B=0, double R_C=0, std::vector<double> P_J={}, bool norm=false);
    void compute_normalisation_factors_excitation_XBC(double E);
    void compute_bound_transition_emission_rates();
    void redistribute_optically_thick_emission(bool uniform = false);
    void match_input_output_rates();
    void compute_continuum_emission_rates();
    void add_line(std::vector<double>& spectrum, double lambda_line, double emission_rate, 
                  double lambda_min, double lambda_max, double d_lambda);
    void add_doppler_broadened_line(std::vector<double>& spectrum, double lambda_line, double emission_rate, 
                                    double T, double lambda_min, double lambda_max, double d_lambda); 




};
#endif // H2SPECTRUM_H