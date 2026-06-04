/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 /* 
This file defines the Precip class, which simulates auroral precipitation using 
a Monte Carlo approach (for the primaries) and two-stream approximation (for the secondaries). 
The class is initialized with source parameters, a source object that defines the initial particle 
distribution, and a reference to the world model for densities and other parameters. 
It loads the necessary cross section data, initializes the simulation variables, 
and provides methods for running the simulation and processing the results.
*/


 // Precip.hpp

#pragma once
#include <vector>
#include <memory>
#include "constants.h"
#include "simparams.h"
#include "Source.hpp"
#include "World1D.hpp"
#include "Egrid.hpp"

class Precip {

    
    Params sp; // source parameters
    std::shared_ptr<Source> src;

    // simulation parameters
    std::vector<float> Zinx;                                // altitude grid (precomputed)
    std::vector<float> Z_edges;                             // altitude bin edges (nz+1) computed from Zinx
    std::vector<float> dzcm;                                 // altitude bin size in cm (precomputed)
    std::vector<float> nH2;                                 // neutral H2 density [cm^-3]
    std::vector<float> sigmasE;                             // cross section for each collision type and energy bin [cm^2]
    std::vector<float> total_sigma;                         // total cross section for each energy bin [cm^2]
    std::vector<float> sigmas;                              // all cross sections for a single energy bin, to be populated in the kernel [cm^2]
    std::vector<float> scatter_angle_pdf;                   // scattering angle PDF for each collision type and energy bin (flattened 2D array: e * nCollisions + c) [unitless]
    std::vector<float> prob;                                // Walker's alias method probability table (flattened 2D array: e * nCollisions + c) [unitless] 
    std::vector<int>   alias;                               // Walker's alias method alias table (flattened 2D array: e * nCollisions + c) [collision type index]
    std::vector<float> deltaEs;                             // Energy loss for each collision type (eV)
    std::vector<std::vector<std::vector<double>>> sigi;     // Cross section into each secondary energy bin for each collision type and primary energy bin (3D array: collision type * primary energy bin * secondary energy bin) [cm^2/eV]
    std::vector<std::vector<double>> secprodsig;            // secondary (actually tertiary) production cross section for each secondary energy bin and collision type (2D array: collision type * secondary energy bin) [cm^2/eV]
    std::vector<float> secprobs;                            // secondary (actually tertiary) production probability for each secondary energy bin and collision type (2D array: collision type * secondary energy bin) [unitless]
    std::vector<std::vector<double>> sigmaiEff;             // effective cross section for inelastic collisions that produce secondaries, for each collision type and energy bin (2D array: collision type * energy bin) [cm^2]
    std::vector<double> sigmaiEffTot;                       // total effective cross section for all inelastic collisions that produce secondaries, for each energy bin (1D array: energy bin) [cm^2]
    std::vector<std::vector<double>> sigSink;               // total removal cross section for each collision type and energy bin


    // simulation variables
    std::vector<float> dt;      // time step in seconds
    std::vector<float> z;       // altitude in m
    std::vector<int> zcell;     // altitude cell index
    std::vector<float> zlocal;  // local z coordinate relative to the cell
    std::vector<float> y;       // horizontal position in m (not to be taken literally due to the absence of a magnetic field)
    std::vector<float> vz;      // velocity in z direction in m/s
    std::vector<float> vy;      // velocity in y direction in m/s
    std::vector<float> E;       // energy in eV
    std::vector<int> alive;     // boolean alive flag (0 or 1)

    //simulation output
    std::vector<int> colcount;                  // collision count for each collision type and energy bin (flattened 2D array: e * nCollisions + c)
    std::vector<int> nion ;                     // number of ions produced (flattened 2D array: z * nE + e) (e at ejected energy)
    std::vector<int> theta_sampled;             // sampled scattering angle for testing
    std::vector<std::vector<double>> qion;      // ionisation rate (2D array: energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    std::vector<std::vector<double>> qion1;     // ionisation rate of primaries only
    std::vector<std::vector<double>> phiPlus;   // upward electron flux as a function of energy and altitude (2D array: energy bin * altitude bin) // cm-2 s-1 eV-1
    std::vector<std::vector<double>> phiMinus;  // downward electron flux as a function of energy and altitude (2D array: energy bin * altitude bin) // cm-2 s-1 eV-1
    std::vector<std::vector<double>> qion2;     // secondary ionisation rate as a function of energy and altitude (2D array: energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    std::vector<std::vector<std::vector<double>>> exRates2; // secondary excitation rates for each collision type, energy bin, and altitude bin (3D array: collision type * energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    std::vector<std::vector<std::vector<double>>> exRates; // total excitation rates for each collision type, energy bin, and altitude bin (3D array: collision type * energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1



public:
    SimParams p;   // simulation parameters

    std::vector<float> qz; // ionisation rate per altitude bin
    std::vector<float> qz1; // ionisation rate per altitude bin
    std::vector<std::vector<double>> exRateB;   // excitation rate for singlet B as a function of energy and altitude (2D array: energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    std::vector<std::vector<double>> exRateC;   // excitation rate for singlet C as a function of energy and altitude (2D array: energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    std::vector<std::vector<double>> exRateEF;  // excitation rate for singlet EF as a function of energy and altitude (2D array: energy bin * altitude bin) // ptle-1 cm-1 s-1 eV-1
    
    // Constructor
    Precip(Params params, std::shared_ptr<Source> src_, const World1D &world_);
    ~Precip() = default;
    void initialisePrimaries() {
        // Keep source sampling bounds aligned with the simulation energy grid.
        sp.Emin = p.E0;
        sp.Emax = p.Emax;
        src->init(sp, z, y, vz, vy, E, alive);
        // Rotate velocities to align with the magnetic field dip angle
        float theta = (90.0f - world.Bdipang) * constants::pi / 180.0f; // angle to rotate by in radians (0 = no rotation = vertical field)
        float vz_old, vy_old;
        for (size_t i = 0; i < sp.N; ++i) {
            vz_old = vz[i];
            vy_old = vy[i];
            vz[i] = vz_old * cos(theta) - vy_old * sin(theta);
            vy[i] = vz_old * sin(theta) + vy_old * cos(theta);
        }
        // Set initial cell-local z coordinates for the stepping algorithm
        for (size_t i = 0; i < sp.N; ++i) {
            zcell[i] = p.Ztoz(z[i]);;
            zlocal[i] = z[i] - p.Z_edges[zcell[i]];
        }
    }

    // Main processing functions
    void processPrimaries();
    void writeSampledPDFToFile();
    void nionToQion();
    void writeQzToFile(bool primariesOnly = false);
    void writeFUVExRatesToFile();
    void writeExRatesToFile();
    void writeColcountToFile();
    void writeQionToFile(bool primariesOnly = false);
    void writeThetaToFile();
    void processSecondaries();
    void sumQionOverE();
    void combinePrimarySecondary();
    void writePhiToFile(std::vector<std::vector<double>>& phi, 
                        const std::string& filename);
    void sampleTheta();
    void run();
    void runlite();


    // Helper functions for coordinate and energy discretization


    // Energy-grid helper owned by Precip so conversion functions use p.nbinsE
    Egrid energyGrid;
    inline double eToE(int e) const { return energyGrid.eToE(e); }      // Convert energy bin index to energy in eV
    inline double eToEcg(int e) const { return energyGrid.eToEcg(e); }  // Convert energy bin index to geometric center energy in eV
    inline double eToEca(int e) const { return energyGrid.eToEca(e); }  // Convert energy bin index to arithmetic center energy in eV
    inline double eTodE(int e) const { return energyGrid.eTodE(e); }    // Convert energy bin index to energy bin width in eV
    inline size_t EToe(double E) const { return static_cast<size_t>(energyGrid.EToIndex(E)); } // Convert energy in eV to energy bin index


 
     // file paths
    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    static inline const std::string datadir = jaicroot + "data/precip/";
    std::string outdir = jaicroot + "out/precip/";
    static inline const std::string nH2vz_file = datadir + "nH2vz.bin";  // H2 number density input data file file

private:

    // Reference to the world model (used for densities, etc.)
    const World1D &world; 

    // Build the alias table for sampling collision types using Walker's method
    void buildAliasTable();

    // Compute the effective cross section for inelastic collisions that produce secondaries, 
    // for each collision type and energy bin 
    void computeSigaEcascade(const std::vector<float>& deltaEs);

    // Compute the backscatter probability for screened rutherford scattering from the PDF
    std::vector<double> backscatterProbabilityRutherford();

    // Solve the tridiagonal system for the two-stream approximation using Thomas algorithm
    void thomasAlgorithm(const std::vector<double>& a, 
                      const std::vector<double>& b, 
                      const std::vector<double>& c, 
                      const std::vector<double>& d, 
                      std::vector<double>& x);

};
