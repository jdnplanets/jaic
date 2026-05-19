/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */


 /*
 This file implements the H2Spectrum class, which calculates the H2 emission spectrum based on excitation rates
and transition probabilities. It reads in transition data, computes excitation rates, and generates the emission spectrum,
 which can be convolved with a Gaussian profile to simulate instrumental broadening. 
 */

 // H2Spectrum.cpp

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>

#include "H2Spectrum.hpp"


using namespace std;
namespace fs = filesystem;

// Global constants
const double k_B_JK = 1.380649e-23;  // Boltzmann constant (J/K)
const double k_B_cm = 0.69503476;    // Boltzmann constant in cm^-1/K
const double c = 2.99792458e8;       // Speed of light (m/s)
const double m_H2 = 3.347e-27;       // Mass of H2 in kg
const double B_e = 60.853;           // Rotational constant in cm^-1 from NIST: https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1000
const double D_e = 0.047;            // Centrifugal distortion constant in cm^-1
const double D00 = 36118.069;        // Dissociation energy in cm^-1

// H2Spectrum class definition
// Constructor: Model initialization (read files, set up normalisation factors, etc.)
H2Spectrum::H2Spectrum() {

    // Here we perform initialization that is not specific to a given run
    // For example, initialize g_J and E_J.

    cout << "Initializing H2Spectrum..." << endl;
    cout << "Reading transition probabilities..." << endl;

    // Read data files
    read_transition_data("B_transitions.dat", "B");   // B<-X R and P branches
    read_transition_data("C-_transitions.dat", "C-"); // C<-X Q
    read_transition_data("C+_transitions.dat", "C+"); // C<-X R and P


    // Filter negligible transitions based on threshold
    Athresh = 5e6;
    filter_negligible_transitions(Athresh);

    cout << "Reading Franck-Condon factors..." << endl;

   
    read_franck_condon_factors_ar("frcosoEFXnc.txt", "E+"); // Abgrall and Roueff for EF<-X O and S branches
    read_franck_condon_factors_ar("frcoqEFXnc.txt", "E-"); // Abgrall and Roueff for EF<-X Q branch
    read_franck_condon_factors_fw("H2_EF1-B1_FCF.dat", "EFB"); // Fantz and Wunderlich data for EF->B

    // Get the maximum J level
    max_Jl = 0;
    for (const auto& dataset : transitions) {
        for (const auto& tr : dataset.second) {
            if (tr.vl == 0) {
                if (tr.Jl > max_Jl) {
                    max_Jl = tr.Jl;
                }
            }
        }
    }

    // Get the maximum J level for the B state
    max_J_B = 0;
    for (const auto& tr : transitions["B"]) {
        if (tr.Ju > max_J_B) {
            max_J_B = tr.Ju;
        }
    }

    // Compute the rovibrational term values and the statistical weights
    for (int J = 0; J <= max_Jl; J++) {
        J_levels.push_back(J);
        double E_J_cm = B_e * J * (J + 1) - D_e * J * J * (J + 1) * (J + 1);
        E_J.push_back(E_J_cm);

        double g_J_value = 2 * J + 1;  // Statistical weight para-H2
        if (J % 2 != 0) {
            g_J_value *= 3; // Ortho-H2
        }
        g_J.push_back(g_J_value); // Statistical weight
    }  

    compute_EF_B_transition_energies();


    cout << "Reading continuum data..." << endl;
    string contdir = datadir + "abgrall_continuum/";
    for (const auto& entry : fs::directory_iterator(contdir)) {
        if (entry.is_regular_file()) {
            loadContinuumFile(entry.path().filename().string(), cmap);
        }
    }

    cout << "...done." << endl;

    // Sum the transition probabilities for each level - bound and continuum
    // for normalisation of the emission step
    sum_Aconts();
    sum_Abounds();



    num_points = static_cast<size_t>((lambda_max - lambda_min) / d_lambda) + 1;

    cout << "H2Spectrum initialized." << endl;

}

// generate() computes the line spectrum for a given set of input parameters.
// T, R_B, R_C, and R_E are provided as inputs.
void H2Spectrum::generate(double T, double E, double R_B, double R_C, double R_E, vector<double>& lambda, vector<double>& spec) {


    auto start0 = chrono::high_resolution_clock::now();
    // Clear the data maps and vectors
    level_excitation_rates.clear();
    level_excitation_rates_ef.clear();
    normalisation_factors_ef.clear();
    normalisation_factors_ef_cascade.clear();
    normalisation_factors_bc_upper.clear();
    normalisation_factors_bc_lower.clear();
    emission_rates.clear();
    continuum_spec.clear();
    lambda_out.clear();
    spec_out.clear();

    total_excitation_rate = R_B + R_C + R_E;
    
    // Compute the normalisation factors for the direct excitation rates and EF cascade
    compute_normalisation_factors_excitation_XEF(E);
    compute_normalisation_factors_cascade_EFB();
    compute_normalisation_factors_excitation_XBC(E);

    // Compute the partition function.
    double Z = compute_partition_function(g_J, E_J, T);

        // Compute the population distribution for the ground state
    vector<double> P_J;
    for (size_t j = 0; j < J_levels.size(); j++) {
        P_J.push_back(g_J[j] * exp(-E_J[j] / (k_B_cm * T)) / Z);
    }

    // Compute the excitation rates for the EF state
    excitation_XEF(E, R_E, P_J);

    // Compute the excitation rates for the B state from the EF to B cascade
    cascade_EFB();

    // Compute the excitation rates for the B and C states
    excitation_XBC(E, R_B, R_C, P_J);

    // Compute the emission rates for each transition
    compute_bound_transition_emission_rates();

    // Redistribute optically thick emission rates
    if (redistribute_thick) redistribute_optically_thick_emission();

    // Compute total excitation rates for B and C states
    double total_emission_rate_B = 0.0;
    double total_emission_rate_C = 0.0;
    for (const auto& entry : emission_rates) {
        const auto& [state, vu, Ju, vl, Jl] = entry.first;
        const EmissionData& data = entry.second;
        if (state == "B") {
            total_emission_rate_B += data.emission_rate;
        } else {
            total_emission_rate_C += data.emission_rate;
        }
    }
    
    // Generate the continuum emission
    compute_continuum_emission_rates();

    // Calculate the total continuum emission
    double total_continuum_emission = accumulate(continuum_spec.begin(), continuum_spec.end(), 0.0);
    totH2 = total_emission_rate_B + total_emission_rate_C + total_continuum_emission;
    match_input_output_rates();
    double total_bound_emission = 0.0;
    for (const auto& entry : emission_rates) total_bound_emission += entry.second.emission_rate;
    total_continuum_emission = accumulate(continuum_spec.begin(), continuum_spec.end(), 0.0);
    totH2 = total_bound_emission + total_continuum_emission;

    // Now put everything on a grid and convolve with a Gaussian if required
    // grid_convolve_spectrum();
    grid_spectrum();

    for (size_t i = 0; i < num_points; i++) {
        lambda[i] = lambda_out[i];
        spec[i] = spec_out[i];
    }

}


// -----------------------------------------------------------
// Data reading functions
// -----------------------------------------------------------

// Function to read transition data, merging C- and C+ into C
void H2Spectrum::read_transition_data(const string& filename, string state) {
    ifstream file(datadir+filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        if (line[0] == '#') continue;

        istringstream iss(line);
        int vu, Ju, vl, Jl;
        double Abound, energy, term_B, term_X;

        if (!(iss >> vu >> Ju >> vl >> Jl >> Abound >> energy >> term_B >> term_X)) {
            cerr << "Error reading line in file: " << filename << endl;
            continue;
        }

        char branch;
        if (Ju == Jl + 1)
            branch = 'R';
        else if (Ju == Jl - 1)
            branch = 'P';
        else
            branch = 'Q';

        transitions[state].push_back({vu, Ju, vl, Jl, Abound, energy, term_B, term_X, branch});
    }
    file.close();
}

// Function to read Franck–Condon factors from Abgrall and Roueff file for EF<-X transitions
void H2Spectrum::read_franck_condon_factors_ar(const string& filename, string state) {
    ifstream file(datadir+filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        if (line[0] == '#') continue;

        istringstream iss(line);
        int vu, Ju, vl, Jl;
        double fcf, energy;

        if (!(iss >> vu >> Ju >> vl >> Jl >> fcf >> energy)) {
            cerr << "Error reading line in file: " << filename << endl;
            continue;
        }

        char branch;
        if (Ju == Jl + 2)
            branch = 'S';
        else if (Ju == Jl - 2)
            branch = 'O';
        else
            branch = 'Q';

        franck_condon_transitions[state].push_back({vu, Ju, vl, Jl, fcf, energy, branch});
    }
    file.close();
}

// Function to read Franck–Condon factors from Fantz and Wunderlich file
void H2Spectrum::read_franck_condon_factors_fw(const string& filename, string state) {
    ifstream file(datadir+filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    vector<vector<double>> fcf_matrix;

    while (getline(file, line)) {
        if (line[0] == '#' || line.find("v'->v''") != string::npos) continue;

        istringstream iss(line);
        vector<double> row;
        double val;
        int v_prime;
        iss >> v_prime;

        while (iss >> val) {
            row.push_back(val);
        }

        fcf_matrix.push_back(row);
    }
    file.close();

    franck_condon_factors[state] = fcf_matrix;
}

// // Function to call the code that reads the continuum data for each file
// void H2Spectrum::read_continuum_data(const string& filename) {
//     cout << "Reading continuum data..." << endl;

//     ifstream infile(datadir+filename);
//     if (!infile.is_open()) {
//         cerr << "Error: Could not open file " << filename << endl;
//         return;
//     }

//     string line;
//     while (getline(infile, line)) {
//         istringstream iss(line);
//         string state;
//         int vu, Ju;
//         double A_tot, df;

//         // Read the line into variables
//         if (iss >> state >> vu >> Ju >> A_tot >> df) {
//             // Store the A_tot value in the map
//             continuum[state][{vu, Ju}] = df;
//         }
//     }

//     infile.close();
// }


void H2Spectrum::filter_negligible_transitions(double threshold) {
    for (auto& dataset : transitions) {
        string state = dataset.first;
        auto& trans_list = dataset.second;

        trans_list.erase(
            std::remove_if(trans_list.begin(), trans_list.end(),
                           [threshold](const Transition& tr) { return tr.Abound < threshold; }),
            trans_list.end());
    }
}

// -----------------------------------------------------------
// Spectrum computation functions
// -----------------------------------------------------------


// Function to compute Hönl–London factors
double H2Spectrum::honl_london(int J, char branch, string state) {

    if (state[0] == 'C'){  // Delta Lambda = 1, Lambda_lower = 0
        if (branch == 'R') {
        return (double) (J + 2);
        } 
        if (branch == 'Q' && J > 0) {
            return ( (2 * J + 1));
        } 
        if (branch == 'P') {
            return (double)(J - 1);
        }
    }
    else {
        if (branch == 'R') { // Delta Lambda = 0
        return (double)(J + 1);
        } 
        if (branch == 'P') {
            return (double)J;
        }
    }
    return 0.0;  // Return 0 if conditions are not met
}


// Function to compute Abgrall et al. (1999) Rotational line strength factors for g-g transitions
double H2Spectrum::Abgrall_SJJ(int Jl, int Ju) {

    double beta = 0.6;
    double term1 = 0.;
    double factor = 0.;
    if (Ju == Jl) {
        term1 = beta;
        factor = Ju * (Ju + 1) / 
                    static_cast<float>((2 * Ju - 1) * (2 * Ju + 3));
    }

    else if (Jl == Ju + 2) {
        factor = 3 * (Ju + 1) * (Ju + 2) / 
                    static_cast<float>(2 * (2 * Ju + 3) * (2 * Ju + 5));
    } 
    else if (Jl == Ju - 2) {
        factor = 3 * Ju * (Ju - 1) / 
                    static_cast<float>(2 * (2 * Ju - 1) * (2 * Ju - 3));
    }
    else factor = 0.;

    return term1 + (1 - beta) * factor;  // Return 0 if conditions are not met
}

// Function to compute Ajello rotational factors
double H2Spectrum::Ajello_SJJ(int Jl, int Ju) {

    double factor = 0.;

    if (Ju == Jl + 1) {
        factor = Ju;
    } 
    else if (Ju == Jl - 1) {
        factor = Ju + 1;
    }
    else factor = 0.;

    return factor;  // Return 0 if conditions are not met
}


// Function to compute Liu factors
double H2Spectrum::Liu_SJJ(int Jl, int Ju) {

    double factor = 0.;

    if (Ju == Jl + 1) {
        factor = std::sqrt(Ju + 1);
    } 
    else if (Ju == Jl - 1) {
        factor = std::sqrt(Ju);
    }
    else factor = 0.;

    return factor;  // Return 0 if conditions are not met
}

// Function to compute the Liu et al. (1995) "effective" rotational line strength correction factor
double H2Spectrum::Liu_correction_factor(const Transition& tr, string state){

    double cp=1.;
    double APJ1=1.;
    double AQJ=1.;

    if (state == "B"){
        // Get AP(J=1)
        for (const auto& tr2 : transitions[state]) {
            if (tr2.vl == tr.vl && tr2.Jl == 1 && tr2.vu == tr.vu && tr2.Ju == 0) {
                APJ1 = tr2.Abound;
                break;
            }
        }
        if (tr.branch == 'P') // P branch
            cp = (2 * tr.Ju + 1) * tr.Abound / ((tr.Ju + 1) * APJ1);
        else
            cp = (2 * tr.Ju + 1) * tr.Abound / (tr.Ju * APJ1);       
    }
    else{
        // Get AQ(J)
        for (const auto& tr2 : transitions[state]) {
            if (tr2.vl == tr.vl && tr2.Jl == tr.Ju && tr2.vu == tr.vu && tr2.Ju == tr.Ju) {
                AQJ = tr2.Abound;
                break;
            }
        }
        if (tr.branch == 'P')
            cp = (2 * tr.Ju + 1) * tr.Abound / (tr.Ju * AQJ);
        else if (tr.branch == 'R')
            cp = (2 * tr.Ju + 1) * tr.Abound / ((tr.Ju + 1) * AQJ);
        else 
            cp = 1.;
    }

    if (isinf(cp) || isnan(cp)) {
        cerr << "Error: Invalid correction factor for state " << state << " and transition " << tr.vl << ", " << tr.Jl << " -> " << tr.vu << ", " << tr.Ju << endl;
        exit(EXIT_FAILURE);
    }

    return cp;

}


// Partition function
double H2Spectrum::compute_partition_function(const vector<double>& g_J, const vector<double>& E_J, double temp) {
    double Z = 0.0;
    for (size_t i = 0; i < g_J.size(); i++) {
        Z += g_J[i] * exp(-E_J[i] / (k_B_cm * temp));
    }
    return Z;
}

// Function to sum continuum transition probabilities for each upper level (v', J')
void H2Spectrum::sum_Aconts(){

    vector<string> states = {"B", "C+","C-"};
    for (const auto& state : states) {
        for (int vu = 0; vu < 39; ++vu) {
            for (int Ju = 0; Ju < 26; ++Ju) {
                for (int Jl = Ju - 1; Jl <= Ju + 1; ++Jl) {

                    vector<double> ek;
                    vector<double> prob;
                    double ev;
                    string str = string(1, state[0]);

                    int res = getContinuumProbForLevel(cmap, str, vu, Ju, Jl, ev, ek, prob, false); // verbose = false

                    if (res == 1) continue;

                    // Integrate prob over ek
                    double integrated_prob = 0.0;
                    for (size_t k = 1; k < ek.size(); k++) {
                        double de = ek[k] - ek[k - 1];
                        double avg_prob = (prob[k] + prob[k - 1]) / 2.0 * de / 2.19474631e+05;  // Convert from cm^-1 to eV
                        integrated_prob += avg_prob;
                    }
                    level_Acont_sums[state][{vu, Ju}] += integrated_prob;
                }
            }
        }
    }
}

// Function to sum bound transition probabilities for each upper level (v', J')
void H2Spectrum::sum_Abounds(){
    // Sum bound transition probabilities for each upper level (v', J')
    for (const auto& dataset : transitions) {
        string state = dataset.first;

        for (const auto& tr : dataset.second) {
                pair<int, int> upper_level = {tr.vu, tr.Ju};
                double cp = Liu_correction_factor(tr, state);
                level_Abound_sums[state][upper_level] += tr.Abound;// * cp;
        }
    }
}

// Function to compute the excitation rates for the X->EF transitions
// Can be run with norm = true to compute the normalisation factors
void H2Spectrum::excitation_XEF(double E, double R_E, vector<double> P_J, bool norm) {

    for (const auto& dataset : franck_condon_transitions) {
        string state = string(1,dataset.first[0]); // Keep the excitation rate over all the EF state
        string state_full = dataset.first; // But keep the full state name for the level excitation rates

        double totexr = 0.0; // Total excitation rate for the current state
        for (const auto& tr : dataset.second) {
            if (tr.vl == 0) {
                pair<int, int> lower_level = {tr.vl, tr.Jl};
                pair<int, int> upper_level = {tr.vu, tr.Ju};
                if (tr.vu >= 9){
                    continue;
                }

                double x = E / (tr.energy/8065.56); //  Both in eV
                double C0 =  0.50490267;
                double C1 = -0.25500813;
                double C2 =  0.24515133;
                double C3 =  0.10720355;
                double C4 = -1.7236746;
                double C5 =  0.41800000;
                double C8 =  0.20983777;
                double C[5] = {C0, C1, C2, C3, C4};
                double term1 = C0 / (x*x) * (1 - 1 / x);
                double term2 = 0.0;
                for (int k = 1; k <= 4; ++k) {
                    term2 += C[k] * (x - 1) * exp(-k * C8 * x);
                }
                double term3 = 1 - 1 / x; // omitting C7 factor in last term
                double Omega_ij = 1/E * C5 * (term1 + term2 + term3);
                // constants get normalised out so we can just keep the variables

                // Calculate the Rotational factor for X -> EF transition
                double S_JJp = Abgrall_SJJ(tr.Jl, tr.Ju);
                double numerator = Omega_ij * tr.fcf * S_JJp;
 

                if (norm){
                    // Accumulate the normalisation factors for the lower level
                    normalisation_factors_ef[state][lower_level] += numerator;
                    level_excitation_rates_ef[state_full][upper_level] = 0.0; // Initialize the excitation rate for this level
                    continue;
                }
                
                // We shouldn't get past here if we are normalising
                double normalisation_factor = normalisation_factors_ef[state][lower_level];

            // Contribution to the excitation rate
                // Ensure we avoid division by zero
                if (normalisation_factor > 0) {

                    // Compute the normalised excitation rate
                    double excitation_rate =
                        R_E * P_J[tr.Jl] * numerator / normalisation_factor;

                    // Accumulate the excitation rate for the upper level
                    level_excitation_rates_ef[state_full][upper_level] += excitation_rate;
                    totexr += excitation_rate;
                }
            }
        }
    }
}
void H2Spectrum::compute_normalisation_factors_excitation_XEF(double E){
    excitation_XEF(E,0,{},true);
}


void H2Spectrum::compute_EF_B_transition_energies(){

    // calculate the excitation rates for the B state from the cascade

    // Loop over all EF vibrational levels and plausible rotational levels
    for (size_t v_ef = 0; v_ef < franck_condon_factors["EFB"].size(); ++v_ef) {
        for (int J_ef = 0; J_ef <= 25; ++J_ef) { // 25 matches other Ju loops used elsewhere
            // Find EF->X transition energy and the corresponding X lower level (if available)
            double EF_X_energy = 0.0;
            int trvl = -1, trJl = -1;
            bool foundEF = false;
            // try E- then E+
            for (const auto& tr : franck_condon_transitions["E-"]) {
                if (tr.vu == static_cast<int>(v_ef) && tr.Ju == J_ef) {
                    EF_X_energy = tr.energy;
                    trvl = tr.vl;
                    trJl = tr.Jl;
                    foundEF = true;
                    break;
                }
            }
            if (!foundEF) {
                for (const auto& tr : franck_condon_transitions["E+"]) {
                    if (tr.vu == static_cast<int>(v_ef) && tr.Ju == J_ef) {
                        EF_X_energy = tr.energy;
                        trvl = tr.vl;
                        trJl = tr.Jl;
                        foundEF = true;
                        break;
                    }
                }
            }
            if (!foundEF) continue;

            // Find the X term value corresponding to that lower X level (from B transition data)
            double term_X = 0.0;
            bool foundXterm = false;
            for (const auto& tr : transitions["B"]) {
                if (tr.vl == trvl && tr.Jl == trJl) {
                    term_X = tr.terml;
                    foundXterm = true;
                    break;
                }
            }
            if (!foundXterm) continue;

            double termEF = EF_X_energy + term_X;

            // For each possible B vibrational level populated by EF->B (from the FCF table)
            for (size_t v_b = 0; v_b < franck_condon_factors["EFB"][v_ef].size(); ++v_b) {
                // Allowed ΔJ = ±1 for EF->B cascade
                for (int J_b : {J_ef - 1, J_ef + 1}) {
                    if (J_b < 0 || J_b > max_J_B) continue;
                    // Find B-state term for (v_b, J_b)
                    double term_B = 0.0;
                    bool foundB = false;
                    for (const auto& tr : transitions["B"]) {
                        if (tr.vu == static_cast<int>(v_b) && tr.Ju == J_b) {
                            term_B = tr.termu;
                            foundB = true;
                            break;
                        }
                    }
                    if (!foundB) continue;
                    double trE = termEF - term_B; // EF->B transition energy
                    if (trE < 0) continue; // Skip unphysical negative transition energies
                    EF_B_transition_energies["EFB"][make_tuple(static_cast<int>(v_ef), J_ef, static_cast<int>(v_b), J_b)] = trE;
                }
            }
        }
    }
}


void H2Spectrum::cascade_EFB(bool norm){

    // Now calculate the excitation rates for the B state from the cascade
    double totexr = 0.0;
    for (const auto& stexr : level_excitation_rates_ef) {
        string state = stexr.first;
        for (const auto& xr : stexr.second) {
            pair<int, int> upper_level = xr.first; // (v_ef, J_ef)
            int v_ef = upper_level.first;  // Vibrational level in the EF state
            int J_ef = upper_level.second; // Rotational level in the EF state

            double efr = xr.second;  // Excitation rate into (v_ef, J_ef)

    
        // For each (v,J) level in the EF state, compute the normalisation factor
        // by summing over all possible decays to the B state
           
            for (int v_b = 0; v_b < franck_condon_factors["EFB"][v_ef].size(); v_b++) {

                double fcf = franck_condon_factors["EFB"][v_ef][v_b];  // Franck-Condon factor for EF -> B

                // ΔJ = ±1, ±3 transitions from the cascade 
                // "three rotational levels, Ji = 1, 3, and 5, can populate the same BSigma1+ level by indirect excitation rate." - Liu 2003
                for (int J_b : {J_ef - 1,  J_ef + 1}) {
                    if (J_b < 0 || J_b > max_J_B) continue;  // Skip unphysical J values
                    
                    char transition_branch = (J_ef == J_b + 1) ? 'R' : 'P';
                    double S_JJp = honl_london(J_b, transition_branch, "B");  // Hönl–London factor

                    auto it_state = EF_B_transition_energies.find("EFB");
                    auto key = make_tuple(v_ef, J_ef, v_b, J_b);
                    if (it_state == EF_B_transition_energies.end()) continue;
                    auto it_trE = it_state->second.find(key);
                    if (it_trE == it_state->second.end()) continue;
                    double trE = it_trE->second;

                    double numerator = fcf * trE * S_JJp;

                    if (norm){
                        normalisation_factors_ef_cascade[state][{v_ef, J_ef}] += numerator;
                        continue;
                    }
                    
                    // We shouldn't get past here if we are normalising
                    double normalisation_factor = normalisation_factors_ef_cascade[state][{v_ef, J_ef}];
                    if (normalisation_factor == 0) continue; 

                    level_excitation_rates["B"][{v_b, J_b}] += efr * numerator / normalisation_factor;
                    totexr += efr * numerator / normalisation_factor;

                    if (isnan(level_excitation_rates["B"][{v_b, J_b}])) {
                        cerr << "Error: NaN detected in level_excitation_rates for v_b = " << v_b << ", J_b = " << J_b << endl;
                        cout << "efr: " << efr << " numerator: " << fcf << " " << S_JJp  << " normalisation factor: " << normalisation_factors_ef_cascade[state][{v_ef, J_ef}] << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            } // End loop over B state levels for this EF level
        } // End loop over EF state levels
    }
}
void H2Spectrum::compute_normalisation_factors_cascade_EFB(){
    // cout << "Computing normalisation factors for the EF->B cascade..." << endl;
    normalisation_factors_ef_cascade.clear();
    cascade_EFB(true);
}


// Function to compute the excitation rates for the X->B and X->C transitions
void H2Spectrum::excitation_XBC(double E, double R_B, double R_C, vector<double> P_J, bool norm) {
    
    for (const auto& dataset : transitions) {
        string state = string(1,dataset.first[0]); // Keep the excitation rate over all the C state
        string state_full = dataset.first; // But keep the full state name for the level excitation rates

        double totsigma = 0.0; // Total cross-section for the current state
        for (const auto& tr : dataset.second) {
            if (tr.vl == 0) {
                pair<int, int> lower_level = {tr.vl, tr.Jl};
                pair<int, int> upper_level = {tr.vu, tr.Ju};

                double x = E / (tr.energy / 8065.56); // Convert Eij from to cm-1 to eV
                double C0 = -0.01555195;
                double C1 = -0.13491574;
                double C2 = -0.02691103;
                double C3 =  0.32786896;
                double C4 = -0.49744809;
                double C5 = -0.43500000;
                double C6 =  0.43500000;
                double C8 =  0.17762538;
                double C[5] = {C0, C1, C2, C3, C4};
                double term1 = C0 / (x*x) * (1 - 1 / x);
                double term2 = 0.0;
                for (int k = 1; k <= 4; ++k) {
                    term2 += C[k] * (x - 1) * exp(-k * C8 * x);
                }
                double term3 = C5 + C6 / x + log(x); // omitting C7 factor in last term
                double Omega_ij = term1 + term2 + term3;
                double gigj = (2 * tr.Ju + 1) / static_cast<float>(2 * tr.Jl + 1); // g_i/g_j
                // constants get normalised out so we can just keep the variables
                double sigma_ij = tr.Abound * pow(tr.energy,-3)/E * gigj * Omega_ij;

                double numerator = sigma_ij;
                

                if (norm){
                    // Accumulate the normalisation factors for the lower level
                    normalisation_factors_bc_lower[state][lower_level] += numerator;
                    continue;
                }
                // We shouldn't get past here if we are normalising
                totsigma += P_J[tr.Jl] * numerator;
                double normalisation_factor = normalisation_factors_bc_lower[state][lower_level];

            // Contribution to the excitation rate
                // Ensure we avoid division by zero
                if (normalisation_factor > 0) {

                    // Compute the normalised excitation rate
                    double excitation_rate = 
                    ((state == "B" ? R_B : R_C) * P_J[tr.Jl]) * numerator / normalisation_factor;// * cp;

                    // Accumulate the excitation rate for the upper level
                    level_excitation_rates[state_full][upper_level] += excitation_rate;

                }
            }
        }
    }

}
void H2Spectrum::compute_normalisation_factors_excitation_XBC(double E){
    // cout << "Computing normalisation factors for the X->B and X->C direct excitation..." << endl;
    excitation_XBC(E,0,0,{},true);
}


// Function to compute the emission rates for each bound transition
void H2Spectrum::compute_bound_transition_emission_rates(){
    emission_rates.clear();

    double emission_rate = 0;
    // Compute emissions for each transition
    for (const auto& dataset : transitions) {
        string state = dataset.first;
        for (const auto& tr : dataset.second) {
            pair<int, int> upper_level = {tr.vu, tr.Ju};
            double level_Atot = level_Abound_sums[state][upper_level] + level_Acont_sums[state][upper_level];

            emission_rate = level_excitation_rates[state][upper_level] * tr.Abound / level_Atot;

            if (emission_rate < threshold) continue; // Skip negligible emission rates
            
            // Add to the emission rates map
            tuple<string, int, int, int, int> key = make_tuple(state, tr.vu, tr.Ju, tr.vl, tr.Jl);
            emission_rates[key].emission_rate += emission_rate;
            emission_rates[key].energy = tr.energy;
            emission_rates[key].wavelength = 1e8 / tr.energy;  // Convert cm^-1 to Å    
        }
    }
}

void H2Spectrum::redistribute_optically_thick_emission(bool uniform){

// Map to store optically thick emission data
    map<string, map<string, map<pair<int, int>, tuple<int, double, double>>>> thick_lines;

    // For each upper level:
    //  - If the lower level is v'' = 0, set the transition counter to 0 and add the emission rate
    //  - Otherwise, increment the transition counter and add the A-coefficient
    for (const auto& dataset : transitions) {
        string state = dataset.first;
        for (const auto& tr : dataset.second) {
            string branch = (tr.Ju == tr.Jl + 1) ? "PR" : (tr.Ju == tr.Jl - 1) ? "PR" : "Q";
            auto& tl = thick_lines[state][branch][{tr.vu, tr.Ju}];
            auto key = make_tuple(state, tr.vu, tr.Ju, tr.vl, tr.Jl);
            // If the emission line exists (i.e. it met the threshold)
            auto it = emission_rates.find(key);
            if (it != emission_rates.end()) {
                const auto& data = it->second;
                // Check if this is an optically thick line
                if (tr.vl == 0){
                    get<0>(tl) = 0;  // Number of transitions
                    get<1>(tl) += data.emission_rate; 
                    get<2>(tl) = 0; // A-coefficient normalisation factor
                }
                else {
                    get<0>(tl) += 1; // Increment the number of transitions
                    get<2>(tl) += tr.Abound; // Increment the A-coefficient normalisation factor
                }
            }
        }
    }

    // Now redistribute the optically thick emission rates
    // Two options: redistribute uniformly or redistribute based on branching ratios
    for (const auto& dataset : transitions) {
        string state = dataset.first;
        for (const auto& tr : dataset.second) {
            auto key = make_tuple(state, tr.vu, tr.Ju, tr.vl, tr.Jl);
            // If the emission line exists (i.e. it met the threshold)
            string branch = (tr.Ju == tr.Jl + 1) ? "PR" : (tr.Ju == tr.Jl - 1) ? "PR" : "Q";
            auto& tl = thick_lines[state][branch][{tr.vu, tr.Ju}];
            auto it = emission_rates.find(key);
            if (it != emission_rates.end()) {
                auto& data = it->second;
                // Check if this is an optically thick line
                if (tr.vl == 0) {
                    // Set the emission rate of this optically thick line to zero (completely reabsorbed)
                    data.emission_rate = 0;
                }
                else{
                    double new_emission_rate = 0;
                    if (uniform)
                    // Redistribute uniformly over the other transitions from this upper level
                        new_emission_rate = data.emission_rate + get<1>(tl) * 1. / static_cast<double>(get<0>(tl));
                    // Redistribute based on branching ratios
                    else 
                        new_emission_rate = data.emission_rate + get<1>(tl) * tr.Abound / static_cast<double>(get<2>(tl));
                    if (new_emission_rate < threshold) continue; // Skip negligible emission rates
                    data.emission_rate = new_emission_rate;
                }
            }
        }
    }
}

void H2Spectrum::compute_continuum_emission_rates(){
// // Initialize the wavelength grid
    continuum_spec.resize(num_points, 0.0);

    double excitation_rate;
    for (const auto& dataset : level_excitation_rates) {
        string state = dataset.first;
        for (const auto& tr : dataset.second) {
            int vu = tr.first.first;
            int Ju = tr.first.second;
            excitation_rate = tr.second;

            pair<int, int> upper_level = {vu, Ju};
            double level_Atot = level_Abound_sums[state][upper_level] + 
                        level_Acont_sums[state][upper_level];

            for (int j = -1; j <= 1; j++) {
                int Jl = Ju + j;
                if (Jl < 0 || Ju > 25 || vu > 28) continue; // catch many transitions with no continuum data to avoid searching for them
                if (state == "B" && Jl == Ju) continue; // Skip transitions to unphysical J values  

                vector<double> ek;
                vector<double> prob;
                double ev;
                int res = getContinuumProbForLevel(cmap, state, vu, Ju, Jl, ev, ek, prob, false);
                if (res == 1){
                    continue;
                }

                double totcont = 0.0;

                for (size_t idx = 1; idx < ek.size(); idx++) {
                    double e = ek[idx];
                    double de = ek[idx] - ek[idx - 1];
                    // Convert to wavelength in Å
                    double lambda = 1e8 / (ev + 118377.2 - (e + D00));

                    // intensity for this energy point
                    double emission_rate = excitation_rate / level_Atot * prob[idx] * de/219474.631;
                    // Add to the continuum spectrum
                    size_t index = static_cast<size_t>((lambda - lambda_min) / d_lambda);
                    // Check that the index is within the grid bounds
                    if (lambda >= lambda_min && lambda <= lambda_max) {
                        continuum_spec[index] += emission_rate;
                        totcont += emission_rate;
                    }

                }
            }
        }
    }
}


void H2Spectrum::match_input_output_rates() {


    double rate_ratio = total_excitation_rate / totH2;
    // Scale the emission rates to match the input excitation rates
    for (auto& entry : emission_rates) {
        entry.second.emission_rate *= rate_ratio;
    }

    for (auto& entry : continuum_spec) {
        entry *= rate_ratio;
    }

}

//--------------------------------------------------------------------------------
// Functions to generate a spectrum on a grid from the lines
//--------------------------------------------------------------------------------

// Function to add a line to the spectrum
void H2Spectrum::add_line(vector<double>& spectrum, double lambda_line, double emission_rate, 
                double lambda_min, double lambda_max, double d_lambda) {
    // Find the closest wavelength index
    size_t index = static_cast<size_t>((lambda_line - lambda_min) / d_lambda);

    // Check that the index is within the grid bounds
    if (lambda_line >= lambda_min && lambda_line <= lambda_max) {
        spectrum[index] += emission_rate;
    }
}

// Function to compute the Doppler-broadened profile for a line and add it to the spectrum
// This is quite slow and not necessary for cold temperatures
void H2Spectrum::add_doppler_broadened_line(vector<double>& spectrum, double lambda_line, double emission_rate, 
                                double T, double lambda_min, double lambda_max, double d_lambda) {
    double sigma_D = lambda_line * sqrt((2 * k_B_JK * T) / m_H2) / c;  // Doppler width in Å
    // cout << "Doppler width for line at " << lambda_line << " Å: " << sigma_D << " Å" << endl;
    for (size_t i = 0; i < spectrum.size(); i++) {
        double lambda = lambda_min + i * d_lambda;
        double gaussian = exp(-pow(lambda - lambda_line, 2) / (2 * pow(sigma_D, 2))) / (sigma_D * sqrt(2 * M_PI));
        spectrum[i] += emission_rate * gaussian;
    }
}

// Function to convolve the spectrum with a Gaussian kernel
vector<double> H2Spectrum::convolve_with_gaussian(vector<double> spectrum, double d_lambda, double fwhm) {
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
            int j = i + k;
            if (j >= 0 && j < spectrum.size()) {
                convolved_spectrum[i] += spectrum[j] * gaussian_kernel[k + kernel_half_width];
            }
        }
    }

    return convolved_spectrum;
}


void H2Spectrum::grid_convolve_spectrum() {

    vector<double> spectrum(num_points, 0.0);
    vector<double> lambda(num_points, 0.0);
    for (size_t i = 0; i < num_points; i++) {
        lambda[i] = lambda_min + i * d_lambda;
        spectrum[i] += continuum_spec[i];
    }
    // first smooth the continuum with FWHM 2.5 A
    spectrum = convolve_with_gaussian(spectrum, d_lambda, 2.5);


    // Add each emission line to the spectrum using a simple line or Doppler broadening
    for (const auto& entry : emission_rates) {
        const EmissionData& data = entry.second;
        double lambda_line = 1e8 / data.energy;  // Convert cm^-1 to Å
        if (lambda_line >= lambda_min && lambda_line <= lambda_max) {
            // add_doppler_broadened_line(spectrum, lambda_line, data.emission_rate, T, lambda_min, lambda_max, d_lambda);
            add_line(spectrum, lambda_line, data.emission_rate, lambda_min, lambda_max, d_lambda);
        }
    }

    // Convolve the spectrum with the Gaussian kernel of specified FWHM
    spectrum = convolve_with_gaussian(spectrum, d_lambda, fwhm);

    for (size_t i = 0; i < num_points; i++) {
        spec_out.push_back(spectrum[i]);
        lambda_out.push_back(lambda[i]);
    }
}

void H2Spectrum::grid_spectrum() {

    vector<double> spectrum(num_points, 0.0);
    vector<double> lambda(num_points, 0.0);
    for (size_t i = 0; i < num_points; i++) {
        lambda[i] = lambda_min + i * d_lambda;
        spectrum[i] += continuum_spec[i];
    }

    // Add each emission line to the spectrum using a simple line or Doppler broadening
    for (const auto& entry : emission_rates) {
        const EmissionData& data = entry.second;
        double lambda_line = 1e8 / data.energy;  // Convert cm^-1 to Å
        if (lambda_line >= lambda_min && lambda_line <= lambda_max) {
            // add_doppler_broadened_line(spectrum, lambda_line, data.emission_rate, T, lambda_min, lambda_max, d_lambda);
            add_line(spectrum, lambda_line, data.emission_rate, lambda_min, lambda_max, d_lambda);
        }
    }

    for (size_t i = 0; i < num_points; i++) {
        spec_out.push_back(spectrum[i]);
        lambda_out.push_back(lambda[i]);
    }
}


void H2Spectrum::print_colour_ratio() {
    // Calculate the total emission between 1550 and 1620 Å
    double cr1 = 0.0;
    for (size_t i = 0; i < num_points; i++) {
        double lambda = lambda_min + i * d_lambda;
        if (lambda >= 1550.0 && lambda <= 1620.0) {
            cr1 += spec_out[i];
        }
    }
    cout << "Total emission between 1550 and 1620 Å: " << cr1 << endl;

    // Calculate the total emission between 1230 and 1300 Å
    double cr2 = 0.0;
    for (size_t i = 0; i < num_points; i++) {
        double lambda = lambda_min + i * d_lambda;
        if (lambda >= 1230.0 && lambda <= 1300.0) {
            cr2 += spec_out[i];
        }
    }
    cout << "Total emission between 1230 and 1300 Å: " << cr2 << endl;
    cout << "Colour ratio: " << cr1 / cr2 << endl;
}

