/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /*
This file defines the H3pColumn class, which models the vertical column of Jupiter's ionosphere 
to compute the H3+ emission spectrum and related quantities. It uses the H3pSpectrum class to 
compute the spectrum based on input parameters such as temperature and H3+ density profiles.
The class also includes methods for applying non-LTE scaling factors, computing volume emission rates, 
and writing output to files. 
*/

 // H3pColumn.hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "H3pSpectrum.hpp"
#include "Source.hpp"
#include "World1D.hpp"


class H3pColumn{
public:

    Params sp;                   // source parameters
    std::shared_ptr<Source> src; // Pointer to the source of excitation rates

    int nz;                      // Number of zinx points
    std::vector<float> Z;        // Altitude grid in km
    float dz;                    // Altitude step in m
    float dzcm;                  // Altitude step in cm
    size_t nw;                   // Number of wavelength points in the output spectrum
    std::vector<double> lambda;  // Wavelength grid in μm
    double obsang = 0;           // Observation angle in degrees (0 = nadir, 90 = limb)
    bool applyNonLTE = true;     // Whether to apply non-LTE scaling factors

    // File paths
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    std::string outdir = jaicroot + "out/h3pspec/";

    // Constructor: Model initialization (read files, set up normalisation factors, etc.)
    H3pColumn(Params params, std::shared_ptr<Source> src_, const World1D &world_, bool applyNonLTE = true);
    ~H3pColumn() = default;


    std::vector<double> specout;    // Emergent spectrum at the top of the column in W m^-2 sr^-1 μm^-1
    std::vector<double> ver;        // volume emission rate vs altitude [W m^-3 sr^-1]
    std::vector<double> q10ver;     // volume emission rate of Q(1,0-) line vs altitude [W m^-3 sr^-1]
    double totef = 0.0;                 // Total observed energy flux in mW m^-2

    // Member functions
    void computeSpectraVsAltitude();
    void computeEmergentSpectrum();
    void computeVERvsAltitude();
    void computeQ10VERvsAltitude();
    double getTotalRadiance();
    double getF335MRadiance();
    int getPeakzIndex();
    void writeVERtoFile();
    void writeQ10VERtoFile();
    void writeEmergentSpectrumToFile();
    void run() {
        computeSpectraVsAltitude();
        computeVERvsAltitude();
        computeQ10VERvsAltitude();
        computeEmergentSpectrum();
        writeVERtoFile();
        writeQ10VERtoFile();
        writeEmergentSpectrumToFile();
    }


protected:

    H3pSpectrum H3pspec;

    // zinx-dependent values
    std::vector<double> n_H3p;                  // H3+ density profile in m^-3
    std::vector<double> T   ;                   // Temperature profile in K    
    std::vector<std::vector<double>> H3pspecvz; // H3+ emission spectrum vs wavelength and altitude, in W m^-2 sr^-1 μm^-1

    const World1D &world;                       // Reference to the world model, for access to Z, T, nH3p, etc.

    // Miller et al. (2013) non-LTE scaling factors
    double getNonLTEScalingFactor(size_t z);
    std::vector<double> nonLTEs;

    // Non-LTE lookup tables (static). Temperatures in K, H2 densities in cm^-3.
    static inline const std::vector<double> nonLTE_temps = {
        300, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000
    };

    static inline const std::vector<double> nonLTE_h2dens = {
        1e12, 1e14, 1e16, 1e18, 1e20
    };

    // Rows correspond to temperatures (same order as nonLTE_temps),
    // columns correspond to H2 densities (same order as nonLTE_h2dens).
    static inline const std::vector<std::vector<double>> nonLTE_scalings = {
        {0.0067, 0.3931, 0.9848, 0.9998, 1.0000},
        {0.0011, 0.0313, 0.7449, 0.9965, 1.0000},
        {0.0013, 0.0108, 0.4955, 0.9866, 0.9999},
        {0.0013, 0.0064, 0.3640, 0.9701, 0.9997},
        {0.0011, 0.0049, 0.2894, 0.9499, 0.9995},
        {0.0010, 0.0042, 0.2469, 0.9299, 0.9992},
        {0.0010, 0.0042, 0.2469, 0.9299, 0.9992},
        {0.0008, 0.0036, 0.2064, 0.9014, 0.9988},
        {0.0007, 0.0034, 0.1961, 0.8923, 0.9987},
        {0.0007, 0.0033, 0.1889, 0.8854, 0.9985},
        {0.0007, 0.0033, 0.1836, 0.8802, 0.9985}
    };


};
