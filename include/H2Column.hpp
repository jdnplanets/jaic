/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines the H2Column class, which computes the H2 emission spectrum from a column of Jupiter's 
atmosphere given excitation rates and atmospheric composition. It includes methods to read excitation rates, 
compute volume emission rates, and calculate the emergent spectrum accounting for absorption by hydrocarbons. 
The class relies on the World1D model for atmospheric properties and a Source for excitation rates. 
*/

 // H2Column.hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "H2Spectrum.hpp"
#include "Source.hpp"
#include "World1D.hpp"
#include "Egrid.hpp"
#include "simparams.h"

class H2Column{
public:

    Params sp; // source parameters
    std::shared_ptr<Source> src; // Pointer to the source of excitation rates

    int nz;                    // Number of zinx points
    std::vector<float> Z;      // Altitude grid in km
    float dz;                  // Altitude step in m
    float dzcm;                // Altitude step in cm
    int nE = 500;               // Number of energy bins for the spectrum calculation (not necessarily the same as the number of energy bins in the excitation rates, which are defined by the source)
    size_t nw;                 // Number of wavelength points in the output spectrum
    std::vector<double> lambda;// Wavelength grid in Å
    bool absorb;               // Flag to include hydrocarbons in the cross sections
    double obsang = 0;         // Observation angle in degrees (0 = zenith, 90 = limb)

    // Output file paths
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    std::string outdir = jaicroot + "out/h2spec/";

    // Constructor: Model initialization (read files, set up normalisation factors, etc.)
    H2Column(Params params, std::shared_ptr<Source> src_, const World1D &world_, const bool hydrocarbon_absorption=true);
    ~H2Column() = default;

    void readExcitationRates(const Params &params, const std::shared_ptr<Source> &src, double F);

    // Set excitation rates for the B, C, and E states.
    // F is an optional flux factor to scale the rates.
    void setExcitationRates(std::vector<std::vector<double>>& R_B,
                            std::vector<std::vector<double>>& R_C, 
                            std::vector<std::vector<double>>& R_E, 
                            const double F= 1.0);


    std::vector<double> specout;    // Emergent spectrum [photons cm-2 s-1 radiating into 4π sr]
    std::vector<double> ver;        // volume emission rate in vs altitude [photons cm-3 s-1]
    std::vector<double> ver_eV;     // volume emission rate vs altitude in eV cm-3 s-1
    double totkr = 0.0;             // Total observed intensity in kRayleighs
    double totkr0 = 0.0;            // Total unabsorbed intensity in kRayleighs observed at zenith
    double cr = 0.0;                // Colour ratio


    // Main processing functions
    void computeSpectraVsAltitude();
    void computeEmergentSpectrum();
    void computeVERvsAltitude();
    void computeVEReVvsAltitude();
    double getTotalUnabsorbedIntensity();
    double getTotalObservedIntensity();
    double getColourRatio();
    void writeVERtoFile();
    void writeVEReVtoFile();
    void writeEmergentSpectrumToFile();
    void run() {
        computeSpectraVsAltitude();
        computeVERvsAltitude();
        computeVEReVvsAltitude();
        computeEmergentSpectrum();
        writeVERtoFile();
        writeVEReVtoFile();
        writeEmergentSpectrumToFile();
    }


protected:

    H2Spectrum h2spec;

    // zinx-dependent values
    std::vector<double> n_H2;               // Neutral H2 density in cm^-3      
    std::vector<double> n_CH4;              // Methane density in cm^-3
    std::vector<double> n_C2H2;             // Acetylene density in cm^-3
    std::vector<double> n_C2H4;             // Ethylene density in cm^-3
    std::vector<double> n_C2H6;             // Ethane density in cm^-3
    std::vector<double> T   ;               // Temperature in K
    std::vector<std::vector<double>> R_B_in ;  // Excitation rates for the B state, in cm^-3 s^-1, indexed by [z][E] read in
    std::vector<std::vector<double>> R_C_in ;  // Excitation rates for the C state, in cm^-3 s^-1, indexed by [z][E] read in
    std::vector<std::vector<double>> R_E_in ;  // Excitation rates for the E state, in cm^-3 s^-1, indexed by [z][E] read in
    std::vector<std::vector<double>> R_B ;  // Excitation rates for the B state, in cm^-3 s^-1, indexed by [z][E]
    std::vector<std::vector<double>> R_C ;  // Excitation rates for the C state, in cm^-3 s^-1, indexed by [z][E]
    std::vector<std::vector<double>> R_E ;  // Excitation rates for the E state, in cm^-3 s^-1, indexed by [z][E]
    std::vector<std::vector<double>> h2specvz; // H2 spectrum vs altitude and wavelength
    std::vector<double> col_CH4;            // Column density of CH4 above each altitude bin in cm^-2
    std::vector<double> col_C2H2;           // Column density of C2H2 above each altitude bin in cm^-2
    std::vector<double> col_C2H4;           // Column density of C2H4 above each altitude bin in cm^-2
    std::vector<double> col_C2H6;           // Column density of C2H6 above each altitude bin in cm^-2

    const World1D &world;
    Egrid energyGrid; // energy grid for excitation bins
    Egrid h2specEGrid; // energy grid for spectrum calculation

    // Structure to hold cross section data for hydrocarbons
    struct CnHnXSectData {
        double lambda;
        double CH4;
        double C2H2;
        double C2H4;
        double C2H6;
    };

    // Vector to hold the hydrocarbon cross section data read from the CSV file
    std::vector<CnHnXSectData> CnHn_xs_data;

    // Function to read hydrocarbon photon cross section data from CSV file
    // Cross sections are in units of cm^2/molecule
    std::vector<CnHnXSectData> readCnHnXSectData(const std::string& filename);

  // A helper structure to hold cross sections for all species.
    struct HydrocarbonCrossSections {
        double CH4;
        double C2H2;
        double C2H4;
        double C2H6;
    };

    // Given a wavelength (in Å) and a vector of cross‐section data (assumed sorted
    // in ascending order by lambda), return the cross sections (for all species) 
    // from the grid point nearest to the requested wavelength.
    HydrocarbonCrossSections getCrossSections(double wavelength,
                                                const std::vector<CnHnXSectData>& xsData)
    ;

    void rebinExcitationRates();


};
