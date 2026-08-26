/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/* 
This file implements the Precip class, which simulates auroral precipitation using 
a Monte Carlo approach (for the primaries) and two-stream approximation (for the secondaries). 
The class is initialized with source parameters, a source object that defines the initial particle 
distribution, and a reference to the world model for densities and other parameters. 
It loads the necessary cross section data, initializes the simulation variables, 
and provides methods for running the simulation and processing the results.
*/



 // Precip.cpp

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <numeric>  // For accumulate
#include <random>
#include <sstream>
#include <iomanip> // For std::setw
#include <algorithm> 
#include <filesystem>
#include <vector>
#include <cmath>
#include <iostream>
#include "Precip.hpp"
#include "FileIO.hpp"
#include "ArrayOps.hpp"
#include "device_arrays.h"


#ifdef USE_CUDA
    #include <cuda_runtime.h> // for CUDA runtime API
#endif

using namespace std;
using namespace constants;


Precip::Precip(Params params, std::shared_ptr<Source> src_, const World1D &world_):
                sp(params), src(std::move(src_)), world(world_)
{


    cout << endl;
    cout << "*************** Starting the auroral precipitation simulation ***************" << endl;
    cout << "Running with the following parameters:" << endl;
    cout << endl;
    cout << "Simulation run ID: " << sp.runid << endl;
    cout << src->banner(sp) << endl;
    cout << "Source label: " << src->label(sp) << endl;
    cout << "Number of primary electrons: " << sp.N << endl;
    cout << "Initial altitude of primary electrons: " << sp.Zinit/1e3 << " km" << endl;
    cout << "Initial energy of primary electrons: " << sp.Einit << " eV" << endl;
    cout << "Total energy in primary electrons: " << sp.N * sp.Einit / 1e3 << " keV" << endl;
    cout << "*****************************************************************************" << endl;
    cout << endl;

    // Populate the simulation parameters with the relevant world and source parameters
    p.P0 =      world.P0;   // bottom of domain in Pa
    p.P1 =      world.P1;   // top of domain in Pa
    p.nbinsP =  world.nP;   // number of pressure bins
    p.Z1 =      world.Z1;   // Top of the simulation domain in m
    p.Z0 =      world.Z0;   // Bottom of the simulation domain in m
    p.nbinsZ =  world.nz;   // Number of altitude bins
    p.N =       sp.N;       // Number of particles
    p.Einit =   sp.Einit;   // Initial energy in eV
    p.Zinit =   sp.Zinit;   // Initial position in m
    p.Z_edges = world.Z_edges.data(); // pointer to altitude bin edges in m, size nz+1
    


    // Initialize energy grid from simulation parameters defined in simparams.h
    p.update_derived();
    energyGrid = Egrid::fromEminEmax(p.nbinsE,static_cast<double>(p.E0),static_cast<double>(p.Emax),p.egridType);
    cout << "Energy grid: " << p.nbinsE << " bins from " << p.E0 << " eV to " 
        << p.Emax << " eV (" << energyGrid.typestring() << ")" << endl;

    // output file path
    outdir += sp.runid+"/";
    // Check if the directory exists, and create it if it doesn't
    if (!std::filesystem::exists(outdir)) {
        if (!std::filesystem::create_directories(outdir)) {
            throw std::runtime_error("Failed to create output directory: " + outdir);
        }
    }

    // Initialize the simulation variables
    dt = vector<float>(sp.N, 1e-10f);
    colcount = vector<int>(p.nCollisions * p.nbinsZ * p.nbinsE, 0); // collision count [c][z][e] for each collision type, energy bin, and altitude bin
    zerr = vector<float>(p.N, 0.0f); // compensation term for z accumulation
    nion = vector<int>(p.nbinsZ * p.nbinsE, 0);
    qion = vector<vector<double>>(p.nbinsE, vector<double>(p.nbinsZ, 0.0));
    qion1 = vector<vector<double>>(p.nbinsE, vector<double>(p.nbinsZ, 0.0));
    qz = vector<float>(p.nbinsZ, 0.0f); // ionisation rate per altitude bin
    qz1 = vector<float>(p.nbinsZ, 0.0f); // ionisation rate from primaries per altitude bin
    exRates = vector<vector<vector<double>>>(p.nCollisions, vector<vector<double>>(p.nbinsZ, vector<double>(p.nbinsE, 0.0))); // collision rates for each collision type, energy bin, and altitude bin (3D array: collision type * energy bin * altitude bin) // s-1 eV-1
    
    // These are for input to other modules
    exRateB = vector<vector<double>>(p.nbinsZ, vector<double>(p.nbinsE, 0.0)); // Excitation rate for singlet B
    exRateC = vector<vector<double>>(p.nbinsZ, vector<double>(p.nbinsE, 0.0)); // Excitation rate for singlet C
    exRateEF = vector<vector<double>>(p.nbinsZ, vector<double>(p.nbinsE, 0.0)); // Excitation rate for singlet E

    
    nH2 = world.nH2;
    Zinx = world.Z;
    dzcm = world.dzcm;
    Z_edges = world.Z_edges;

    // Define vector that holds the cross section data for all collisions
    string sigfile;
    if (p.decades == 7.0f) {
        sigfile = "sigmas_all_BEB_SLARex_Rel_10MeV.bin";
    } else {
        sigfile = "sigmas_all_BEB_SLARex_Rel.bin";
    }
    cout << "Loading cross section data from " << sigfile << endl;
    string sigmasE_file = datadir + sigfile;
    sigmasE = vector<float>(p.nE * p.nCollisions, 0.0f);  // cross section data - all cross sections for all runtime energies
    std::vector<float> sigmasE_in = std::vector<float>(p.nE_in * p.nCollisions, 0.0f);  // cross sections on input grid
    utils::io::readFloatBinaryfileData(sigmasE_file, 
                            &sigmasE_in[0], 
                            p.nE_in, 
                            p.nCollisions);

    
    // Saved data are in 10^-16 cm^2. convert to cm^2
    for (int i = 0; i < p.nE_in; i++) {
        for (int j = 0; j < p.nCollisions; j++) {
            sigmasE_in[i * p.nCollisions + j] *= 1e-16; // cm^2
        }
    }

    sigmas = vector<float>(p.nCollisions, 0.0f);  // cross section data - all cross sections for a single energy to be populated in the kernel    


    // Now, if necessary, we need to interpolate the cross section data from the input energy grid 
    // (which is logarithmic with 6 decades and 100 bins) onto the runtime energy grid 
    // (which may be linear or logarithmic, depending on the parameters).

    // If the runtime energy grid is linear, resample the cross sections
    // (which are defined on a 6-decade logarithmic grid) onto the linear grid.
    if (energyGrid.type == EgridType::Linear) {
        Egrid EGrid_in(p.nE_in, p.decades_in, p.E0_in, p.egridType_in);
        std::vector<float> sig_interp(sigmasE.size(), 0.0f);
        for (int e = 0; e < p.nE; ++e) {
            // energy at linear bin center
            double E_lin = energyGrid.eToEcg(e);
            int idx = EGrid_in.EToIndex(E_lin);
            int idx2 = std::min(idx + 1, p.nE_in - 1);
            double Ea = EGrid_in.eToE(idx);
            double Eb = EGrid_in.eToE(idx2);
            double w = 0.0;
            if (Eb > Ea) w = (E_lin - Ea) / (Eb - Ea);
            for (int c = 0; c < p.nCollisions; ++c) {
                double sa = sigmasE_in[idx * p.nCollisions + c];
                double sb = sigmasE_in[idx2 * p.nCollisions + c];
                double val = (1.0 - w) * sa + w * sb;
                sig_interp[e * p.nCollisions + c] = static_cast<float>(val);
            }
        }
        sigmasE.swap(sig_interp);

    } else {
        // Non-linear energy grid: interpolate sigmasE from input grid to runtime grid
        if (p.nE != p.nE_in) {
            Egrid EGrid_in(p.nE_in, p.decades_in, p.E0_in, p.egridType_in);
            std::vector<float> sig_interp(sigmasE.size(), 0.0f);
            for (int e = 0; e < p.nE; ++e) {
                // energy at runtime bin center
                double E = energyGrid.eToEcg(e);
                int idx = EGrid_in.EToIndex(E);
                int idx2 = std::min(idx + 1, p.nE_in - 1);
                double Ea = EGrid_in.eToE(idx);
                double Eb = EGrid_in.eToE(idx2);
                double w = 0.0;
                if (Eb > Ea) w = (E - Ea) / (Eb - Ea);
                for (int c = 0; c < p.nCollisions; ++c) {
                    double sa = sigmasE_in[idx * p.nCollisions + c];
                    double sb = sigmasE_in[idx2 * p.nCollisions + c];
                    double val = (1.0 - w) * sa + w * sb;
                    sig_interp[e * p.nCollisions + c] = static_cast<float>(val);
                }
            }
            sigmasE.swap(sig_interp);
        } else {
            // Runtime energy grid is the same as input energy grid: just copy the data
            sigmasE = sigmasE_in;
        }

    }

    // Calculate the total cross section for each energy
    total_sigma = vector<float>(p.nE, 0.0f);
    for (int i = 0; i < p.nE; i++) {
        float sigsum = 0.0f;
        for (int j = 0; j < p.nCollisions; j++) {
        
            sigsum += sigmasE[i * p.nCollisions + j];
        }
        total_sigma[i] = sigsum;
    }


    string scatter_file;
    if (p.decades == 7.0f) {
        scatter_file = "scatter_angle_pdf_1k_10MeV.bin";
    } else {
        scatter_file = "scatter_angle_pdf_1k.bin";
    } 
    cout << "Loading scatter angle data from " << scatter_file << endl;
    string scatter_angle_pdf_file = datadir + scatter_file;
    scatter_angle_pdf = vector<float>(p.nTheta * p.nE, 0.0f);
    std::vector<float> scatter_angle_pdf_in = std::vector<float>(p.nTheta * p.nE_in, 0.0f);
    utils::io::readFloatBinaryfileData(scatter_angle_pdf_file, 
                            &scatter_angle_pdf_in[0], 
                            p.nTheta, 
                            p.nE_in);
    
                  
    // Resample scatter-angle PDF onto linear energy grid if required
    if (energyGrid.type == EgridType::Linear) {
        Egrid EGrid_in(p.nE, 6.0, 1.0, EgridType::Logarithmic);
        std::vector<float> scat_interp(scatter_angle_pdf.size(), 0.0f);
        for (int t = 0; t < p.nTheta; ++t) {
            for (int e = 0; e < p.nE; ++e) {
                double E_lin = energyGrid.eToEcg(e);
                int idx = EGrid_in.EToIndex(E_lin);
                int idx2 = std::min(idx + 1, p.nE - 1);
                double Ea = EGrid_in.eToE(idx);
                double Eb = EGrid_in.eToE(idx2);
                double w = 0.0;
                if (Eb > Ea) w = (E_lin - Ea) / (Eb - Ea);
                double sa = scatter_angle_pdf_in[t * p.nE_in + idx];
                double sb = scatter_angle_pdf_in[t * p.nE_in + idx2];
                double val = (1.0 - w) * sa + w * sb;
                scat_interp[t * p.nE + e] = static_cast<float>(val);
            }
        }
        scatter_angle_pdf.swap(scat_interp);
    } else {
        // Non-linear energy grid: interpolate sigmasE from input grid to runtime grid
        if (p.nE != p.nE_in) {
            Egrid EGrid_in(p.nE_in, p.decades_in, p.E0_in, p.egridType_in);
            std::vector<float> scat_interp(scatter_angle_pdf.size(), 0.0f);
            for (int t = 0; t < p.nTheta; ++t) {
                for (int e = 0; e < p.nE; ++e) {
                    double E = energyGrid.eToEcg(e);
                    int idx = EGrid_in.EToIndex(E);
                    int idx2 = std::min(idx + 1, p.nE_in - 1);
                    double Ea = EGrid_in.eToE(idx);
                    double Eb = EGrid_in.eToE(idx2);
                    double w = 0.0;
                    if (Eb > Ea) w = (E - Ea) / (Eb - Ea);
                    double sa = scatter_angle_pdf_in[t * p.nE_in + idx];
                    double sb = scatter_angle_pdf_in[t * p.nE_in + idx2];
                    double val = (1.0 - w) * sa + w * sb;
                    scat_interp[t * p.nE + e] = static_cast<float>(val);
                }
            }
            scatter_angle_pdf.swap(scat_interp);
        } else {
            // Runtime energy grid is the same as input energy grid: just copy the data
            scatter_angle_pdf = scatter_angle_pdf_in;
        }
    }



    theta_sampled = vector<int>(p.nTheta * p.nE, 0.0f);
    buildAliasTable();
    // sampleTheta();
    // writeThetaToFile();


    deltaEs = vector<float>(p.nCollisions, 0.0f); //eV

    deltaEs[p.elastic]                = 0.0f;
    deltaEs[p.ionisation]             = 15.43f; // ionisation potential in eV
    deltaEs[p.excitation_triplet_a]   = 11.8f;
    deltaEs[p.excitation_triplet_b]   = 13.2f;
    deltaEs[p.excitation_triplet_c]   = 11.9f;
    deltaEs[p.excitation_triplet_e]   = 13.9f;
    deltaEs[p.excitation_singlet_B]   = 11.2f;
    deltaEs[p.excitation_singlet_C]   = 12.3f;
    deltaEs[p.vibrational]            = 0.516f;
    deltaEs[p.rotational]             = 0.04f;
    deltaEs[p.excitation_singlet_EF]  = 12.6f;


    initialisePrimaries();
    writeSampledPDFToFile();

    // Check that the maximum electron energy is less than the Egrid max
    float maxE = *std::max_element(E.begin(), E.end());
    if (maxE > p.Emax) {
        cout << "WARNING: Maximum electron energy (" << maxE << " eV) exceeds energy grid maximum (" 
             << p.Emax << " eV)" << endl;
    }

    // Initialise redistribution matrix for secondary electron energy cascade
    computeSigaEcascade(deltaEs);   

} // End of Precip constructor


// Build Walker's alias tables for a 2D probability distribution.
// Inputs:
//   scatter_angle_pdf - Flattened vector (size nTheta * nE) where each row corresponds to Theta and
//                       each column corresponds to an energy E.
//   p.nTheta            - Number of Theta (vertical) bins (rows).
//   p.nE                - Number of energy (horizontal) bins (columns).
// Outputs:
//   prob  - Flattened vector (size nTheta * nE) of adjusted probabilities for each (Theta, E).
//   alias - Flattened vector (size nTheta * nE) of alias indices for each (Theta, E).
void Precip::buildAliasTable() {
    // Resize output arrays.

    prob.resize(p.nTheta * p.nE);
    alias.resize(p.nTheta * p.nE);

    // Process each column (each energy value).
    // Since the input is stored in row-major order,
    // the element at row t and column e is at index: t * p.nE + e.
    for (int e = 0; e < p.nE; ++e) {
        // Extract the column corresponding to energy e.
        vector<float> colWeights(p.nTheta);
        float colSum = 0.0f;
        for (int t = 0; t < p.nTheta; ++t) {
            float w = scatter_angle_pdf[t * p.nE + e];
            colWeights[t] = w;
            colSum += w;
        }
        if (colSum <= 0.0f) {
            throw runtime_error("Sum of probabilities in a column is zero or negative.");
        }

        // Normalise the weights so that their sum equals p.nTheta.
        vector<float> norm(p.nTheta);
        for (int t = 0; t < p.nTheta; ++t) {
            norm[t] = colWeights[t] * p.nTheta / colSum;
        }

        // Partition indices into two lists: small (< 1.0) and large (>= 1.0).
        vector<int> small, large;
        for (int t = 0; t < p.nTheta; ++t) {
            if (norm[t] < 1.0f)
                small.push_back(t);
            else
                large.push_back(t);
        }

        // Build the alias table for this column.
        while (!small.empty() && !large.empty()) {
            int l = small.back();
            small.pop_back();
            int g = large.back();
            large.pop_back();

            // For the (Theta, E) position corresponding to row l:
            // Save the adjusted probability and alias index.
            prob[l * p.nE + e] = norm[l];
            alias[l * p.nE + e] = g;

            // Adjust the large index's normalized probability.
            norm[g] = norm[g] - (1.0f - norm[l]);
            if (norm[g] < 1.0f)
                small.push_back(g);
            else
                large.push_back(g);
        }

        // For any indices remaining (either in large or small), set probability to 1.
        while (!large.empty()) {
            int g = large.back();
            large.pop_back();
            prob[g * p.nE + e] = 1.0f;
        }
        while (!small.empty()) {
            int l = small.back();
            small.pop_back();
            prob[l * p.nE + e] = 1.0f;
        }
    }
}


// Compute the histogram of velocities sampled by the source and write to a file
void Precip::writeSampledPDFToFile() {

    vector<int> pdfE = vector<int>(p.nE, 0);
    for (int i = 0; i < p.N; ++i) {
        int k = p.EToe(E[i]);
        pdfE[k] += 1;
    }
    stringstream ss;
    string suf = src->label(sp);
    ss << outdir << "pdf_" << suf << ".dat";
    string pdffilename = ss.str();
    ofstream pdffile(pdffilename, ios::trunc);
    if (pdffile.is_open()) {   
            std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "sampled PDF(E)"},
        {"units", " eV-1"},
        {"source_label", suf},
        {"nbinsE", std::to_string(p.nbinsE)}
    };

    const std::vector<std::string> cols = {"E [eV]", "PDF(E) [eV^-1]"};

    const bool ok = utils::io::write_dat_table_fixed_width(
        pdffilename,
        "Precip model input: E",
        meta,
        cols,
        p.nbinsE,
        [&](int e, std::ostream& os, int w) {
            os << std::right << std::setw(w) << p.El(e)
            << " "        << std::setw(w) << pdfE[e];
        },
        /*col_width=*/12,
        /*precision=*/6
    );

    if (ok) std::cout << "Wrote sampled PDF to " << pdffilename << "\n";
    }

}


void runPrimariesKernel(DeviceArrays darrs, DeviceArrays harrs, SimParams p);

void Precip::processPrimaries(){
    

    cout << "Processing primaries..." << endl;

    DeviceArrays harrs;
    harrs.nH2 = &nH2[0];
    harrs.sigmasE = &sigmasE[0];
    harrs.total_sigma = &total_sigma[0];
    harrs.sigmas = &sigmas[0];
    harrs.prob = &prob[0];
    harrs.alias = &alias[0];
    harrs.dt = &dt[0];
    harrs.z = &z[0];
    harrs.zerr = &zerr[0];
    harrs.y = &y[0];
    harrs.vz = &vz[0];
    harrs.vy = &vy[0];
    harrs.E = &E[0];
    harrs.alive = &alive[0];
    harrs.nion = &nion[0];
    harrs.theta_sampled = &theta_sampled[0];
    harrs.colcount = &colcount[0];


#ifdef USE_CUDA
    DeviceArrays darrs;
    // Allocate and copy host data to device for each array

    size_t size_Zinx   = Zinx.size() * sizeof(float);
    cudaMalloc((void**)&darrs.Zinx, size_Zinx);
    cudaMemcpy(darrs.Zinx, Zinx.data(), size_Zinx, cudaMemcpyHostToDevice);

    size_t size_nH2    = nH2.size() * sizeof(float);
    cudaMalloc((void**)&darrs.nH2, size_nH2);
    cudaMemcpy(darrs.nH2, nH2.data(), size_nH2, cudaMemcpyHostToDevice);

    size_t size_sigmasE    = sigmasE.size() * sizeof(float);
    cudaMalloc((void**)&darrs.sigmasE, size_sigmasE);
    cudaMemcpy(darrs.sigmasE, sigmasE.data(), size_sigmasE, cudaMemcpyHostToDevice);

    size_t size_totalSigma = total_sigma.size() * sizeof(float);
    cudaMalloc((void**)&darrs.total_sigma, size_totalSigma);
    cudaMemcpy(darrs.total_sigma, total_sigma.data(), size_totalSigma, cudaMemcpyHostToDevice);

    size_t size_sigmas     = sigmas.size() * sizeof(float);
    cudaMalloc((void**)&darrs.sigmas, size_sigmas);
    // sigmas will be updated in the kernel if needed; copy if required:
    cudaMemcpy(darrs.sigmas, sigmas.data(), size_sigmas, cudaMemcpyHostToDevice);

    size_t size_prob       = prob.size() * sizeof(float);
    cudaMalloc((void**)&darrs.prob, size_prob);
    cudaMemcpy(darrs.prob, prob.data(), size_prob, cudaMemcpyHostToDevice);

    size_t size_alias      = alias.size() * sizeof(int);
    cudaMalloc((void**)&darrs.alias, size_alias);
    cudaMemcpy(darrs.alias, alias.data(), size_alias, cudaMemcpyHostToDevice);

    size_t size_deltaEs    = deltaEs.size() * sizeof(float);
    cudaMalloc((void**)&darrs.deltaEs, size_deltaEs);
    cudaMemcpy(darrs.deltaEs, deltaEs.data(), size_deltaEs, cudaMemcpyHostToDevice);


    size_t size_dt         = dt.size() * sizeof(float);
    cudaMalloc((void**)&darrs.dt, size_dt);
    cudaMemcpy(darrs.dt, dt.data(), size_dt, cudaMemcpyHostToDevice);

    // Result arrays

    size_t size_z          = z.size() * sizeof(float);
    cudaMalloc((void**)&darrs.z, size_z);
    cudaMemcpy(darrs.z, z.data(), size_z, cudaMemcpyHostToDevice);
    
    size_t size_zerr       = zerr.size() * sizeof(float);
    cudaMalloc((void**)&darrs.zerr, size_zerr);
    cudaMemcpy(darrs.zerr, zerr.data(), size_zerr, cudaMemcpyHostToDevice);

    size_t size_y          = y.size() * sizeof(float);
    cudaMalloc((void**)&darrs.y, size_y);
    cudaMemcpy(darrs.y, y.data(), size_y, cudaMemcpyHostToDevice);

    size_t size_vz         = vz.size() * sizeof(float);
    cudaMalloc((void**)&darrs.vz, size_vz);
    cudaMemcpy(darrs.vz, vz.data(), size_vz, cudaMemcpyHostToDevice);

    size_t size_vy         = vy.size() * sizeof(float);
    cudaMalloc((void**)&darrs.vy, size_vy);
    cudaMemcpy(darrs.vy, vy.data(), size_vy, cudaMemcpyHostToDevice);

    size_t size_E          = E.size() * sizeof(float);
    cudaMalloc((void**)&darrs.E, size_E);
    cudaMemcpy(darrs.E, E.data(), size_E, cudaMemcpyHostToDevice);

    size_t size_alive      = alive.size() * sizeof(int);
    cudaMalloc((void**)&darrs.alive, size_alive);
    cudaMemcpy(darrs.alive, alive.data(), size_alive, cudaMemcpyHostToDevice);

    size_t size_nion       = nion.size() * sizeof(int);
    cudaMalloc((void**)&darrs.nion, size_nion);
    cudaMemcpy(darrs.nion, nion.data(), size_nion, cudaMemcpyHostToDevice);

    size_t size_theta     = theta_sampled.size() * sizeof(int);
    cudaMalloc((void**)&darrs.theta_sampled, size_theta);
    cudaMemcpy(darrs.theta_sampled, theta_sampled.data(), size_theta, cudaMemcpyHostToDevice);

    size_t size_colcount   = colcount.size() * sizeof(int);
    cudaMalloc((void**)&darrs.colcount, size_colcount);
    cudaMemcpy(darrs.colcount, colcount.data(), size_colcount, cudaMemcpyHostToDevice);

    // Copy simulation parameter arrays to device
    float* d_Z_edges = nullptr;
    size_t size_Z_edges   = world.Z_edges.size() * sizeof(float);
    cudaMalloc(&d_Z_edges, size_Z_edges);
    cudaMemcpy(d_Z_edges, world.Z_edges.data(), size_Z_edges, cudaMemcpyHostToDevice);

    SimParams dp = p; // Make a copy of the SimParams struct to pass to the kernel
    dp.Z_edges = d_Z_edges; // Update the pointer in the struct to point to the device memory

    // Launch the CUDA kernel
    runPrimariesKernel(darrs, harrs, dp);

    cudaFree(d_Z_edges); // Needs to be here because the symbol is here
#endif
}


void Precip::nionToQion() {
    // Convert nion in number bin-1 cm-2 to qion in cm-3 eV-1
    // A cm-2 column from the MC run is assumed
    // nion  = nion(E->, z)
    // Transpose nion to [e][z] because in the solver below we pass over z for each energy bin
    for (int e = 0; e < p.nbinsE; e++) {
        double dE = p.dE(e);
        for (int z = 0; z < p.nbinsZ; z++) {
            int thisnion = nion[z * p.nbinsE + e];
            qion1[e][z] = thisnion / static_cast<float>(sp.N) / dzcm[z] / dE; // ptle-1 cm-1 s-1 eV-1
        }
    }

    double max_qz = 0.0;
    for (int z = 0; z < p.nbinsZ; ++z) {
        double sum = 0.0;
        for (int e = 0; e < p.nbinsE; ++e) {
            sum += qion[e][z];
        }
        if (sum > max_qz) max_qz = sum;
    }
    // std::cout << "Maximum value of qz*1e2 from primaries: " << max_qz*1e2 << std::endl;
}


void Precip::writeColcountToFile() {
    string colcount_file = outdir + "colcount.dat";
    ofstream outfile(colcount_file, ios::trunc);
    if (outfile.is_open()) {
        for (int i = 0; i < colcount.size(); ++i) {
            outfile << colcount[i] << " ";
            if ((i + 1) % p.nbinsE == 0)
                outfile << endl;
        }
        outfile.close();
        cout << "Colcount data written to " << colcount_file << endl;
    } else {
        cout << "Unable to open file for writing colcount data" << std::endl;
    }
}


#include "Rnd.h"
Rnd rnd; // Global random number generator
void Precip::sampleTheta(){

    // For every energy bin
    for (int e = 0; e < p.nE; ++e) {


        // Sample theta using the Walker alias method
        // Generate random integers index j in the range [0, nTheta)
        int nRands = 1000000;
        for (int i = 0; i < nRands; i++) {
            int j = static_cast<int>(rnd() * p.nTheta);
            // Generate a random number in [0,1]
            float r = rnd();
            // Get the theta value
            float theta;
            int t;
            if (r < prob[j * p.nE + e]) {
                // theta = (1 - j / static_cast<float>(p.nTheta)) * pi;
                t = j;
            } else {
                // theta = (1 - alias[j * p.nE + e] / static_cast<float>(p.nTheta)) * pi;
                t = alias[j * p.nE + e];
            }
            theta_sampled[(p.nTheta - 1 - t) * p.nE + e]++;
        }
    }

}


void Precip::writeThetaToFile() {

    // Normalize each column of theta_sampled by the total number of counts in each energy bin
    for (int e = 0; e < p.nE; ++e) {
        int total_counts = 0;
        for (int t = 0; t < p.nTheta; ++t) {
            total_counts += theta_sampled[t * p.nE + e];
        }
        if (total_counts > 0) {
            for (int t = 0; t < p.nTheta; ++t) {
                theta_sampled[t * p.nE + e] /= (static_cast<float>(total_counts)/1000000.);
            }
        }
    }
    ofstream thetafile(outdir + "theta_sampled.dat", ios::trunc);
    if (thetafile.is_open()) {
        for (int i = 0; i < p.nTheta; i++) {
            for (int j = 0; j < p.nE; j++) {
                thetafile << theta_sampled[i * p.nE + j] << " ";
            }
            thetafile << endl;
        }
        thetafile.close();
    } else {
        cout << "Unable to open file for writing theta data" << endl;
    }
    cout << "Theta sampling data written to " << outdir + "theta_sampled.dat" << endl;
}



//**************************************************************************************
//Secondary electron processes
//**************************************************************************************
// Energy-cascade probability maps
// - Fills:
//     sigi[c][j][k]    : per-event probability primary goes from bin i to k for channel c
//     secprodsig[j][et]: per-event probability secondary (from ionization) goes into et
//     sigSink[c][j] : per-event probability of "E < T" events for channel c and initial bin i
//**************************************************************************************

void Precip::computeSigaEcascade(const std::vector<float>& deltaEs)
{
    // ----- Energy grid: edges, widths, centers -----
    std::vector<double> Eedge(p.nE + 1), El(p.nE), Eh(p.nE), dE(p.nE), Ec(p.nE);
    for (int j = 0; j <= p.nE; ++j) Eedge[j] = p.El(j);
    for (int j = 0; j < p.nE; ++j) {
        El[j] = Eedge[j];
        Eh[j] = Eedge[j + 1];
        dE[j] = Eh[j] - El[j];
        Ec[j] = p.Eca(j);
    }

    // sigi[collision][initial_bin][target_bin]
    sigi.assign(p.nCollisions,
                std::vector<std::vector<double>>(p.nE, std::vector<double>(p.nE, 0.0)));
    secprodsig.assign(p.nE, std::vector<double>(p.nE, 0.0));

    // sigSink[collision][initial_bin] : effective cross section for "E < T" events
    sigSink.assign(p.nCollisions, std::vector<double>(p.nE, 0.0));


    auto overlap_len = [](double a0, double a1, double b0, double b1) -> double {
        const double lo = std::max(a0, b0);
        const double hi = std::min(a1, b1);
        return std::max(0.0, hi - lo);
    };

    // ----- Build sigi for inelastic (non-ionisation) discrete losses -----
    for (int c = 0; c < p.nCollisions; ++c) {

        double T = static_cast<double>(deltaEs[c]);

        if (c == p.elastic || c == p.ionisation) continue;

        
        // ---- Loop over initial energy bins JY (0-based: j) ----
        for (int i = 0; i < p.nE; ++i) {

            const double sig = sigmasE[i * p.nCollisions + c];
            if (sig <= 0.0) continue;

            // --- Split bin i into below-threshold and above-threshold parts ---
            // Below-threshold: E in [El[i], min(Eh[i], T)) => Sink event (remove from cascade)
            // Above-threshold: E in [max(El[i], T), Eh[i]] => normal fixed-loss T cascade

            const double El_i = El[i];
            const double Eh_i = Eh[i];

            // Entire bin below threshold -> all such "collisions" are Sink events.
            if (Eh_i <= T) {
                sigSink[c][i] += sig;                     // removal cross section for this channel+bin
                continue;
            }

            // Partial below-threshold portion
            double below_len = std::max(0.0, std::min(Eh_i, T) - El_i);
            if (below_len > 0.0) {
                const double frac_below = below_len / dE[i];
                const double sig_below  = sig * frac_below;

                sigSink[c][i] += sig_below;
            }

            // Above-threshold part handled by Swartz+overlap
            const double El_eff = std::max(El_i, T);
            const double Eh_eff = Eh_i;
            const double above_len = Eh_eff - El_eff;
            if (above_len <= 0.0) continue;

            const double frac_above = above_len / dE[i];
            const double sig_above  = sig * frac_above;

            // Effective initial energy for the participating (above-threshold) sub-interval.
            // Use linear mean because you're doing overlap in linear energy.
            const double Ei_eff = 0.5 * (El_eff + Eh_eff);

            // Degraded interval from the above-threshold sub-interval after losing T
            double a = El_eff - T;
            double b = Eh_eff - T;

            // Clamp to >= 0 for indexing/overlap
            if (a < 0.0) a = 0.0;
            if (b <= a) continue;


            // If the whole interval was below El[0], dump all to sink
            if (b <= El[0]) {
                sigSink[c][i] += sig_above;  // remove those events entirely
                continue;
            }

            // Candidate target-bin indices from the (in-grid) degraded interval [a,b]
            size_t kl = p.EToe(a);
            size_t kh = p.EToe(b);

            // Enforce k < i (Swartz requirement Ek < Ei implies strictly lower bin)
            size_t im1 = static_cast<size_t>(i - 1);
            if (kh > im1) kh = im1;
            if (kl > im1) kl = im1;

            // Save original degraded interval before any clamping to El[0]
            const double a0 = a;
            const double b0 = b;

            // ---- siphon sub-grid part (< El[0]) into Sink sink instead of bin 0 ----
            below_len = 0.0;
            if (a0 < El[0]) {
                below_len = std::max(0.0, std::min(b0, El[0]) - a0);
            }

            // Total span (should be > 0 here)
            const double span = b0 - a0;
            if (span <= 0.0) continue;

            // Fraction of events whose post-loss energy would be below the grid
            const double f_below = below_len / span;

            // Add that fraction to the sink (terminal removal)
            if (f_below > 0.0) {
                sigSink[c][i] += sig_above * f_below;
            }

            // Now restrict the interval to the in-grid part for overlap weights
            double a_in = std::max(a0, El[0]);
            double b_in = b0;
            if (b_in <= a_in) {
                // Everything went below grid, nothing to redistribute
                continue;
            }

            // Candidate target-bin indices from the in-grid degraded interval [a_in, b_in]
            kl = p.EToe(a_in);
            kh = p.EToe(b_in);

            // Enforce k < i (Swartz requirement)
            im1 = static_cast<size_t>(i - 1);
            if (kh > im1) kh = im1;
            if (kl > im1) kl = im1;

            // Build overlap weights across [kl..kh] using ONLY in-grid overlaps
            std::vector<double> wk;
            wk.reserve(kh - kl + 1);

            double wsum = 0.0;
            for (size_t k = kl; k <= kh; ++k) {
                const double ol = overlap_len(a_in, b_in, El[k], Eh[k]);
                wk.push_back(ol);
                wsum += ol;
            }

            // If overlaps vanished, dump remaining (in-grid) fraction to nearest lower bin
            if (wsum <= 0.0) {
                const size_t k = im1;
                const double Ek = Ec[k];
                const double denom = Ei_eff - Ek;
                if (denom <= 0.0) continue;

                // Only the remaining fraction (1 - f_below) is redistributed
                const double sigmaEff = (sig_above * (1.0 - f_below)) * (T / denom);
                sigi[c][i][k] += sigmaEff * (dE[i] / dE[k]);
                continue;
            }

            // IMPORTANT: The redistribution only applies to the remaining fraction (1 - f_below)
            const double frac_in = (1.0 - f_below);

            // Normalize weights and compute Ek as weighted mean of the group
            for (double &w : wk) w /= wsum;

            double Ek = 0.0;
            for (size_t idx = 0; idx < wk.size(); ++idx) {
                Ek += wk[idx] * Ec[kl + idx];
            }

            // Enforce Ek < Ei (strictly)
            if (!(Ek < Ei_eff)) Ek = Ec[i - 1];      // or clamp to something < Ei_eff
            const double denom = Ei_eff - Ek;
            if (denom <= 0.0) continue;
            const double sigmaEff = sig_above * frac_in * (T / denom);


            // Spread over target group with dE ratio
            for (size_t idx = 0; idx < wk.size(); ++idx) {
                const size_t k = kl + idx;
                sigi[c][i][k] += sigmaEff * wk[idx] * (dE[i] / dE[k]);
            }
        }
    }

// ######################################################################

    // -------------------------
    // Ionisation contribution
    // -------------------------
    const bool do_ionisation = true;

        if (do_ionisation) {
        const double B = 8.3;

        auto mass_int = [&](double x) -> double {
            if (x <= 0.0) return 0.0;
            return B * std::atan(x / B);
        };
        auto first_moment_int = [&](double x) -> double {
            if (x <= 0.0) return 0.0;
            const double u = x / B;
            return 0.5 * B * B * std::log(1.0 + u * u);
        };

        const int cIon = p.ionisation;
        const double ionT = static_cast<double>(deltaEs[cIon]); // ionisation threshold/potential [eV]
        if (ionT <= 0.0) return; // no ionisation, nothing to do

        // Below-grid secondary sink diagnostics (eps < El[0]) for each primary bin i
        std::vector<double> secIon_sig_below(p.nE, 0.0);  // [cm^2]
        std::vector<double> secIon_eps_below(p.nE, 0.0);  // [cm^2 * eV]

        // Loop over initial energy bins i
        for (int i = 0; i < p.nE; ++i) {

            const double sigIonTot = static_cast<double>(sigmasE[i * p.nCollisions + cIon]);
            if (sigIonTot <= 0.0) continue;

            const double Eci  = Ec[i];
            const double dEi  = dE[i];

            // Fortran-like Tmax = (Eci - ionT)/2, but do not proceed if non-positive
            double Tmax = 0.5 * (Eci - ionT);
            if (Tmax <= 0.0) continue;
            if (Tmax > 1.0e6) Tmax = 1.0e6; // optional cap (kept from your earlier version)

            // Probability normalization over [0, Tmax]
            const double Z = mass_int(Tmax);
            if (Z <= 0.0) continue;

            const double EminGrid = El[0]; // energy-grid minimum

            // below-grid secondaries eps in [0, min(EminGrid, Tmax)) go to sink 
            const double Ecut = std::min(EminGrid, Tmax);
            if (Ecut > 0.0) {
                const double m_lo  = mass_int(Ecut) - mass_int(0.0);
                const double fm_lo = first_moment_int(Ecut) - first_moment_int(0.0);

                if (m_lo > 0.0) {
                    const double prob_lo = m_lo / Z;
                    const double sig_lo  = sigIonTot * prob_lo;   // [cm^2]
                    const double eps_mean_lo = fm_lo / m_lo;      // [eV]

                    secIon_sig_below[i] += sig_lo;
                    secIon_eps_below[i] += sig_lo * eps_mean_lo;
                }
            }

            // Uppermost representable secondary bin index
            int itmax = p.EToe(Tmax);
            if (itmax < 0) continue;
            if (itmax >= p.nE) itmax = p.nE - 1;

            // Loop over secondary-energy target bins k (representable only)
            for (int k = 0; k <= itmax; ++k) {

                // Secondary energy segment [E1, E2] within [EminGrid, Tmax]
                double E1 = std::max(El[k], EminGrid);
                double E2 = Eh[k];
                if (E2 > Tmax) E2 = Tmax;
                if (E2 <= E1) continue;

                // Probability mass in this secondary segment
                const double m12 = mass_int(E2) - mass_int(E1);
                if (m12 <= 0.0) continue;

                const double prob = m12 / Z;
                const double sigion = sigIonTot * prob; // [cm^2] partial ionisation XS into this segment

                // Mean ejected energy in this segment
                const double fm12 = first_moment_int(E2) - first_moment_int(E1);
                double eps_mean = fm12 / m12;
                if (!(eps_mean >= 0.0)) eps_mean = 0.0;

                // --- Secondary production operator (SEC analogue)
                secprodsig[i][k] += sigion * (dEi / dE[k]);

                // ------------------------------------------------------------
                // Primary energy-loss redistribution for this partial channel
                // split sub-threshold to sink, siphon below-grid degraded to sink,
                // Swartz+overlap on the remainder)
                // ------------------------------------------------------------
                const double dE1 = ionT + eps_mean; // energy removed from primary in this partial channel
                if (dE1 <= 0.0) continue;

                const double El_j = El[i];
                const double Eh_j = Eh[i];

                // Entire primary bin below threshold -> all these ionisation events Sink the primary
                if (Eh_j <= dE1) {
                    sigSink[cIon][i] += sigion;
                    continue;
                }

                // Partial below-threshold portion in the primary bin
                const double len_below_thr = std::max(0.0, std::min(Eh_j, dE1) - El_j);
                if (len_below_thr > 0.0) {
                    const double frac_below_thr = len_below_thr / dE[i];
                    sigSink[cIon][i] += sigion * frac_below_thr;
                }

                // Above-threshold portion that can undergo full loss dE1
                const double El_eff = std::max(El_j, dE1);
                const double Eh_eff = Eh_j;
                const double len_above_thr = Eh_eff - El_eff;
                if (len_above_thr <= 0.0) continue;

                const double frac_above_thr = len_above_thr / dE[i];
                const double sigion_above = sigion * frac_above_thr;

                // Effective initial energy for participating portion of primary bin
                const double Ei_eff = 0.5 * (El_eff + Eh_eff);

                // Degraded interval from the above-threshold part after losing dE1
                double a = El_eff - dE1;
                double b = Eh_eff - dE1;

                if (a < 0.0) a = 0.0;
                if (b <= a) continue;

                // ---- siphon the part of [a,b] that lies below El[0] into the sink ----
                const double a0 = a;
                const double b0 = b;
                const double span = b0 - a0;
                if (span <= 0.0) continue;

                const double len_below_grid =
                    (a0 < El[0]) ? std::max(0.0, std::min(b0, El[0]) - a0) : 0.0;

                const double f_below_grid = len_below_grid / span;

                if (f_below_grid > 0.0) {
                    sigSink[cIon][i] += sigion_above * f_below_grid;
                }

                const double frac_in = 1.0 - f_below_grid;
                if (frac_in <= 0.0) continue;

                // Restrict overlap interval to the on-grid part
                const double a_in = std::max(a0, El[0]);
                const double b_in = b0;
                if (b_in <= a_in) continue;

                // Candidate target-bin indices from the in-grid degraded interval [a_in, b_in]
                size_t kl = p.EToe(a_in);
                size_t kh = p.EToe(b_in);

                // Enforce target strictly below i
                const size_t im1 = static_cast<size_t>(i - 1);
                if (kh > im1) kh = im1;
                if (kl > im1) kl = im1;

                // Build overlap weights across [kl..kh] using ONLY in-grid overlaps
                std::vector<double> wk;
                wk.reserve(kh - kl + 1);

                double wsum = 0.0;
                for (size_t kk = kl; kk <= kh; ++kk) {
                    const double Lk = overlap_len(a_in, b_in, El[kk], Eh[kk]);
                    wk.push_back(Lk);
                    wsum += Lk;
                }
                

                // If overlaps vanished, dump remaining fraction to nearest lower bin
                if (wsum <= 0.0) {
                    if (i == 0) {
                        sigSink[cIon][i] += sigion_above * frac_in;
                        continue;
                    }

                    const size_t kfb = im1;
                    const double Ek_fallback = Ec[kfb];
                    const double denom = Ei_eff - Ek_fallback;
                    if (denom <= 0.0) continue;

                    const double sigmaEff = (sigion_above * frac_in) * (dE1 / denom);
                    sigi[cIon][i][kfb] += sigmaEff * (dE[i] / dE[kfb]);
                    continue;
                }

                for (double &w : wk) w /= wsum;

                double Ek = 0.0;
                for (size_t idx = 0; idx < wk.size(); ++idx) {
                    Ek += wk[idx] * Ec[kl + idx];
                }

                if (!(Ek < Ei_eff)) Ek = Ec[i - 1];

                const double denom = Ei_eff - Ek;
                if (denom <= 0.0) continue;

                const double sigmaEff = (sigion_above * frac_in) * (dE1 / denom);

                for (size_t j = 0; j < wk.size(); ++j) {
                    const size_t kk = kl + j;
                    sigi[cIon][i][kk] += sigmaEff * wk[j] * (dE[i] / dE[kk]);
                }
            }
        }
    }



    // ----- Precompute effective loss sigmas for L -----
    sigmaiEff.assign(p.nE, std::vector<double>(p.nCollisions, 0.0));
    sigmaiEffTot.assign(p.nE, 0.0);

    for (int c = 0; c < p.nCollisions; ++c) {
        if (c == p.elastic) continue;

        for (int e = 0; e < p.nE; ++e) {
            double sum = 0.0;

            if (e == 0) {
                sigmaiEff[0][c] = sigi[c][0][0];
            } else {
                for (int k = 0; k < e; ++k) sum += sigi[c][e][k] * dE[k];
                sigmaiEff[e][c] = sum / dE[e];
            }

            sigmaiEff[e][c] += sigSink[c][e];
            sigmaiEffTot[e] += sigmaiEff[e][c];
        }
    }

}



vector<double> Precip::backscatterProbabilityRutherford(){
    
    vector<double> cumsump(p.nTheta * p.nE, 0.0f);
    for (int e = 0; e < p.nE; ++e) {
        double sum = 0.0f;
        for (int t = 0; t < p.nTheta; ++t) {
            sum += scatter_angle_pdf[t * p.nE + e];
            cumsump[t * p.nE + e] = sum;
        }
        // Normalize the cumulative sum for this energy bin
        for (int t = 0; t < p.nTheta; ++t) {
            cumsump[t * p.nE + e] /= sum;
        }
    }

    vector<double> backscatterSums(p.nE, 0.0f);
    for (int e = 0; e < p.nE; ++e) {
        backscatterSums[e] = cumsump[(p.nTheta / 2 - 1) * p.nE + e];
    }

    return backscatterSums;

}

// Helper to compute derivatives
// Uses central difference apart from at the edges, where forward/backward difference is used
void compute_derivative(const std::vector<double>& arr, const std::vector<float>& dz, std::vector<double>& deriv) {
    size_t n = arr.size();
    for (size_t z = 0; z < n; ++z) {
        if (z == 0) {
            deriv[z] = (arr[z + 1] - arr[z]) / dz[z];
        } else if (z == n - 1) {
            deriv[z] = (arr[z] - arr[z - 1]) / dz[z - 1];
        } else {
            deriv[z] = (arr[z + 1] - arr[z - 1]) / (2.0 * dz[z]);
        }
    }
}


// Helper function to refine the vertical grid by a factor of zfact

struct ZGrid {
    std::vector<float> Z;        // lower edges, size nbinsZ
    std::vector<float> Z_edges;  // all edges, size nbinsZ + 1
    std::vector<float> dz;       // cell widths in m, size nbinsZ
    std::vector<float> dzcm_z;    // vertical cell widths in cm, size nbinsZ
    std::vector<float> dzcm_s;    // path lengthcell widths in cm, size nbinsZ
};

inline ZGrid refineZ(
    const std::vector<float>& Z_edges_in,
    int zfact,
    float sinpsi
) {
    const int nz_coarse = static_cast<int>(Z_edges_in.size()) - 1;

    if (nz_coarse <= 0) {
        throw std::runtime_error("refineZ: need at least one cell");
    }

    if (zfact <= 0) {
        throw std::runtime_error("refineZ: zfact must be positive");
    }

    ZGrid out;

    const int nz_fine = nz_coarse * zfact;

    out.Z.reserve(static_cast<std::size_t>(nz_fine));
    out.Z_edges.reserve(static_cast<std::size_t>(nz_fine + 1));
    out.dz.reserve(static_cast<std::size_t>(nz_fine));
    out.dzcm_z.reserve(static_cast<std::size_t>(nz_fine));
    out.dzcm_s.reserve(static_cast<std::size_t>(nz_fine));

    for (int z = 0; z < nz_coarse; ++z) {
        const float z0 = Z_edges_in[z];
        const float z1 = Z_edges_in[z + 1];

        const float dz_coarse = z1 - z0;

        if (!(dz_coarse > 0.0f)) {
            throw std::runtime_error(
                "refineZ: input Z_edges must be strictly increasing"
            );
        }

        const float dz_fine = dz_coarse / static_cast<float>(zfact);

        for (int j = 0; j < zfact; ++j) {
            const float z_lower = z0 + static_cast<float>(j) * dz_fine;

            out.Z.push_back(z_lower);
            out.Z_edges.push_back(z_lower);
            out.dz.push_back(dz_fine);
            out.dzcm_z.push_back(dz_fine * 1.0e2f);
            out.dzcm_s.push_back(dz_fine * 1.0e2f / sinpsi);
        }
    }

    // Final right-hand edge of the top bin
    out.Z_edges.push_back(Z_edges_in.back());

    return out;
}


// ***************************************************************************************
// ***************************************************************************************

void Precip::processSecondaries() {

    cout << "Processing secondaries..." << endl;


    const double psi_rad = world.Bdipang * M_PI / 180.0;
    const double sinpsi = std::sin(psi_rad);

    if (sinpsi <= 0.0) {
        throw std::runtime_error(
            "processSecondaries: sin(psi) must be > 0. "
            "psi = 90 deg corresponds to vertical transport."
        );
    }

    double total_qion_E = 0.0;  
    vector<double> qion_E(p.nbinsE, 0.0); // qion integrated over vertical z for each energy bin

    for (int e = 0; e < p.nbinsE; ++e) {
        double dE = p.dE(e);
        double Ec = p.Ecg(e);
        double E = p.El(e);

        for (int z = 0; z < p.nbinsZ; ++z) {
            total_qion_E += qion1[e][z] * dE * dzcm[z] * Ec;
            qion_E[e]    += qion1[e][z] * dzcm[z];
        }
    }

    cout << "Total integrated qion energy over energy and height per input particle: "
         << total_qion_E << " eV" << std::endl;


    // Convert sigmasE to a 2D array for convenience
    vector<vector<double>> sigsE(p.nbinsE, vector<double>(p.nCollisions, 0.0));

    for (int e = 0; e < p.nbinsE; ++e) {
        for (int c = 0; c < p.nCollisions; ++c) {
            sigsE[e][c] = static_cast<double>(sigmasE[e * p.nCollisions + c]);
        }
    }


    // -------------------------------------------------------------------------
    // Refined vertical grid
    // -------------------------------------------------------------------------

    const int zfactint = 5;
    float zfact = static_cast<float>(zfactint);

    ZGrid zg_fine = refineZ(Z_edges, zfactint, sinpsi);

    const int nz = static_cast<int>(zg_fine.Z.size());

    // Lower edges, size nz
    std::vector<float>& Z_fine = zg_fine.Z;

    // Full edges, size nz + 1
    std::vector<float>& Z_edges_fine = zg_fine.Z_edges;

    // Widths, size nz
    std::vector<float>& dz = zg_fine.dz;
    std::vector<float>& dzcm_z = zg_fine.dzcm_z;
    std::vector<float>& dzcm_s = zg_fine.dzcm_s; //Fine path-length spacing along tilted field/beam



    // Interpolate world.nH2 onto the refined vertical z grid
    std::vector<float> nH2 = utils::array::interp<float>(Zinx, world.nH2, Z_fine);

    // Compute the backscatter probability for Rutherford scattering as a function of energy
    vector<double> rutherfordBSP = backscatterProbabilityRutherford();


    // Compute L(E,z) = n(σ_i,total + p_e σ_e)
    // and     S(E,z) = n(p_e σ_e)

    vector<vector<double>> L, S; 
    L.resize(p.nbinsE, vector<double>(nz, 0.0));
    S.resize(p.nbinsE, vector<double>(nz, 0.0));

    for (size_t e = 0; e < p.nbinsE; ++e) {

        float sigma_e = sigsE[e][p.elastic]; // cm2

        double L_n = sigmaiEffTot[e] + rutherfordBSP[e] * sigma_e; // cm2
        double S_n = rutherfordBSP[e] * sigma_e;                   // cm2

        for (size_t z = 0; z < nz; ++z) {
            L[e][z] = L_n * nH2[z]; 
            S[e][z] = S_n * nH2[z];
        }
    }


    double mu = p.cosThetaMean; // mu = <cos(theta)>


    vector<vector<double>> qp;
    vector<vector<double>> qm;

    qp.resize(p.nbinsE, vector<double>(nz, 0.0));
    qm.resize(p.nbinsE, vector<double>(nz, 0.0));


    // Output arrays
    // phiPlus[e][z], phiMinus[e][z] because innermost loop is over z
    phiPlus.resize(p.nbinsE, vector<double>(nz, 0.0));
    phiMinus.resize(p.nbinsE, vector<double>(nz, 0.0));

    qion2.resize(p.nbinsE, vector<double>(nz, 0.0));

    // Secondary excitation arrays use [c][z][E] to match primaries array order
    exRates2.resize(p.nCollisions, vector<vector<double>>(nz, vector<double>(p.nbinsE, 0.0)));

    double ediss = 0.0;   // energy dissipation/sink integrated over path
    double eoutTop = 0.0; // energy out the top, integrated over E

    //########################################################################################

    // Main loop: solve the parabolic equation, cascade, and compute secondary
    // production for each energy bin.
    //
    // Start from the highest energy bin and work downwards.
    // Safe unsigned-downwards loop avoids underflow when e reaches 0.
    for (size_t e = p.nbinsE; e-- > 0; ) {

        vector<double> q(nz, 0.0);

        // These are derivatives with respect to the path coordinate s, not z.
        vector<double> dLds(nz, 0.0);
        vector<double> dSds(nz, 0.0);
        vector<double> dqmds(nz, 0.0);
        vector<double> dqds(nz, 0.0);

        double dE = p.dE(e); // eV

        double sig_B  = sigsE[e][p.excitation_singlet_B];
        double sig_C  = sigsE[e][p.excitation_singlet_C];
        double sig_EF = sigsE[e][p.excitation_singlet_EF];


        // Interpolate q(E) onto refined vertical z grid.
        vector<double> qe(p.nbinsZ, 0.0);

        for (size_t z = 0; z < p.nbinsZ; ++z) {
            qe[z] = qion1[e][z];
        }

        std::vector<float> qf = utils::array::interp<float>(Zinx, qe, Z_fine);
        q.assign(qf.begin(), qf.end());

        // Compute gradients with respect to the tilted path coordinate s.
        //
        // Since neighbouring altitude samples are separated by dz, but the
        // transport coordinate separation is ds = dz/sin(psi), using dzcm_s
        // here gives d/ds = sin(psi) d/dz.
        compute_derivative(L[e],  dzcm_s, dLds);
        compute_derivative(S[e],  dzcm_s, dSds);
        compute_derivative(q,     dzcm_s, dqds);
        compute_derivative(qm[e], dzcm_s, dqmds);


        // Coefficients for the tridiagonal system
        vector<double> a(nz, 0.0);
        vector<double> b(nz, 0.0);
        vector<double> c(nz, 0.0);
        vector<double> d(nz, 0.0);

        double g0 = 0.0;
        double b0 = 0.0;

        for (size_t z = 0; z < nz; ++z) {

            double alpha = -1.0 / S[e][z] * dSds[z];

            double beta =
                -1.0 / mu *
                (
                    dLds[z]
                    + L[e][z] * L[e][z] / mu
                    - S[e][z] * S[e][z] / mu
                    - L[e][z] / S[e][z] * dSds[z]
                );

            double gamma =
                1.0 / mu *
                (
                    S[e][z] * qp[e][z] / mu
                    + L[e][z] * qm[e][z] / mu
                    + q[z] / (2.0 * mu) * (L[e][z] + S[e][z])
                    + alpha * (q[z] / 2.0 + qm[e][z])
                    + dqds[z] / 2.0
                    + dqmds[z]
                );

            if (z == 0) {
                g0 = gamma;
                b0 = beta;
            }

            a[z] = 1.0 + alpha * dzcm_s[z] / 2.0;
            b[z] = beta * dzcm_s[z] * dzcm_s[z] - 2.0;
            c[z] = 1.0 - alpha * dzcm_s[z] / 2.0;
            d[z] = -gamma * dzcm_s[z] * dzcm_s[z];
        }


        // Solve the tridiagonal system for phiMinus using the Thomas algorithm
        vector<double> phim(nz, 0.0);

        phim[0] = -g0 / b0; // Boundary condition at z = 0, derivatives = 0

        thomasAlgorithm(a, b, c, d, phim);


        // Integrate phiPlus along the tilted path coordinate s
        vector<double> phip(nz, 0.0);

        phip[0] = phim[0]; // Bottom boundary condition

        for (size_t z = 0; z < nz - 1; ++z) {

            double A = L[e][z] / mu;
            double B = S[e][z] * phim[z] / mu
                     + q[z] / (2.0 * mu)
                     + qp[e][z] / mu;

            double x = A * dzcm_s[z];

            if (x > 50.0) x = 50.0; // Limit optical depth to avoid overflow

            phip[z + 1] = (phip[z] - B / A) * exp(-x) + B / A;
        }


        // ---------------------------------------------------------------------
        // Energy cascade for this energy bin
        // ---------------------------------------------------------------------

        double edissE = 0.0;
        double totqm = 0.0;
        double totqp = 0.0;

        // Electrons falling into lower bins due to inelastic collisions.
        // For each inelastic collision, for each target energy bin below this one,
        // compute the resulting production rate at all altitudes.
        for (size_t coll = 1; coll < p.collisionCount; ++coll) { // Exclude elastic collisions

            double p_i; // Backscatter probability

            if (coll == p.ionisation ||
                coll == p.excitation_singlet_B || 
                coll == p.excitation_singlet_C ||
                coll == p.excitation_singlet_EF) {
                p_i = rutherfordBSP[e];
            }
            else {
                p_i = 0.5; // isotropic scattering
            }

            for (size_t z = 0; z < nz; ++z) {
                for (size_t et = 0; et < e; ++et) {

                    const double sig = sigi[coll][e][et];

                    if (sig == 0.0) continue;

                    qp[et][z] += nH2[z] *
                                 (
                                     p_i * sig * phim[z]
                                     + (1.0 - p_i) * sig * phip[z]
                                 );

                    qm[et][z] += nH2[z] *
                                 (
                                     p_i * sig * phip[z]
                                     + (1.0 - p_i) * sig * phim[z]
                                 );
                }
            }
        }


        // Path-integrated energy dissipation.
        for (int coll = 0; coll < p.collisionCount; ++coll) {

            if (coll == p.elastic) continue;

            for (size_t z = 0; z < nz; ++z) {

                double delE = deltaEs[coll];

                edissE += nH2[z] * (phim[z] + phip[z]) * dE * dzcm_s[z]
                        * sigsE[e][coll] * delE;

                ediss  += nH2[z] * (phim[z] + phip[z]) * dE * dzcm_s[z]
                        * sigsE[e][coll] * delE;
            }
        }


        // ---------------------------------------------------------------------
        // Secondary ionisation production
        // ---------------------------------------------------------------------

        double totqsec = 0.0;

        for (size_t et = 0; et < e; ++et) {

            const double sec = secprodsig[e][et];

            if (sec == 0.0) continue;

            for (size_t z = 0; z < nz; ++z) {

                const double phi_sum = phim[z] + phip[z];
                const double sq = nH2[z] * phi_sum * sec;

                qion2[et][z] += sq;

                qm[et][z] += 0.5 * sq;
                qp[et][z] += 0.5 * sq;
            }
        }


        // Store phiPlus and phiMinus
        for (size_t z = 0; z < nz; ++z) {
            phiPlus[e][z]  = phip[z];
            phiMinus[e][z] = phim[z];
        }


        // Secondary excitation rates.
        for (size_t coll = 0; coll < p.collisionCount; ++coll) {

            if (coll == p.elastic) continue;

            for (size_t z = 0; z < nz; ++z) {
                exRates2[coll][z][e] = nH2[z] * (phim[z] + phip[z]) * dE * sigsE[e][coll];
            }
        }


        // Additional sink term, path-integrated along path
        for (int coll = 0; coll < p.collisionCount; ++coll) {

            if (coll == p.elastic) continue;

            for (size_t z = 0; z < nz; ++z) {

                ediss += nH2[z] * (phim[z] + phip[z]) * dE * dzcm_s[z]
                       * sigSink[coll][e] * p.Ec(e);
            }
        }


        // Energy escaping through the top.
        eoutTop += phip[nz - 1] * mu * p.Ec(e) * dE;

    } // end energy loop ###########################################################


    double E_low = 0.0;

    for (size_t z = 0; z < nz; ++z) {
        E_low += (phiMinus[0][z] + phiPlus[0][z])
               * p.dE(0) * p.Ec(0) * dzcm_s[z];
    }
    std::cout << "Energy dissipated: " << ediss << " eV" << std::endl;
    std::cout << "Energy out the top: " << eoutTop << " eV" << std::endl;
    std::cout << "Total energy accounted for (dissipation + energy sink + outflow): "
              << (ediss + eoutTop) << " eV" << std::endl;
    

    // -------------------------------------------------------------------------
    // Reinterpolate qion2, phiPlus, phiMinus, exRates2 back onto the original vertical z axis.
    // -------------------------------------------------------------------------
    std::vector<std::vector<double>> qion2_coarse;
    std::vector<std::vector<double>> phiPlus_coarse; 
    std::vector<std::vector<double>> phiMinus_coarse;
    std::vector<std::vector<std::vector<double>>> exRates2_coarse;

    qion2_coarse.resize(p.nbinsE, vector<double>(p.nbinsZ, 0.0));
    phiPlus_coarse.resize(p.nbinsE, vector<double>(p.nbinsZ, 0.0));
    phiMinus_coarse.resize(p.nbinsE, vector<double>(p.nbinsZ, 0.0));
    exRates2_coarse.resize(p.nCollisions, vector<vector<double>>(p.nbinsZ, vector<double>(p.nbinsE, 0.0)));


    for (int e = 0; e < p.nbinsE; ++e) {
        qion2_coarse[e] =
            utils::array::interp<double>(Z_fine, qion2[e], Zinx);
        phiPlus_coarse[e] =
            utils::array::interp<double>(Z_fine, phiPlus[e], Zinx);
        phiMinus_coarse[e] =
            utils::array::interp<double>(Z_fine, phiMinus[e], Zinx);
    }

    std::vector<double> er(nz, 0.0);
    for (int coll = 0; coll < p.collisionCount; ++coll) {
        for (int e = 0; e < p.nbinsE; ++e) {
            for (int z = 0; z < nz; ++z) {
                er[z] = exRates2[coll][z][e];
            }
            std::vector<double> er_coarse =
                utils::array::interp<double>(Z_fine, er, Zinx);
            for (int z = 0; z < p.nbinsZ; ++z) {
                exRates2_coarse[coll][z][e] = er_coarse[z];
            }
        }
    }


    qion2 = std::move(qion2_coarse);
    phiPlus = std::move(phiPlus_coarse);
    phiMinus = std::move(phiMinus_coarse);
    exRates2 = std::move(exRates2_coarse);


} // End of processSecondaries

// Thomas algorithm for solving tridiagonal system
void Precip::thomasAlgorithm(const vector<double>& a, 
                             const vector<double>& b, 
                             const vector<double>& c, 
                             const vector<double>& d, 
                             vector<double>& x) {
    int n = d.size();
    vector<double> l(n, 0.0);
    vector<double> k(n, 0.0);

    // Forward elimination
    k[0] = (d[1] - c[1] * x[0]) / b[1];
    l[0] = a[1] / b[1];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] - c[i] * l[i-1];
        k[i] = (d[i] - c[i] * k[i-1]) / denom;
        l[i] = a[i] / denom;
    }

    // BC sets derivative at the top to zero
    x[n - 2] = k[n - 2];
    x[n - 1] = x[n - 2];

    // Back substitution
    for (int i = n - 3; i >= 0; --i) {
        x[i] = k[i] - l[i] * x[i+1];
    }
}



void Precip::combinePrimarySecondary() {

    // Combine the primary and secondary ionisation rates
    for (int e = 0; e < p.nbinsE; e++) {
        for (int z = 0; z < p.nbinsZ; z++) {
            qion[e][z] = qion1[e][z] + qion2[e][z]; // cm-3 s-1 eV-1
        }
    }

    //Combine primary and secondary excitation rates
    for (int coll = 0; coll < p.nCollisions; ++coll) {
        for (int z = 0; z < p.nbinsZ; ++z) {
            for (int e = 0; e < p.nbinsE; ++e) {
                int thisnex = colcount[coll * p.nbinsZ * p.nbinsE + z * p.nbinsE + e];
                float exRate1 = static_cast<float>(thisnex) / static_cast<float>(sp.N) / dzcm[z]; // cm-3 s-1
                exRates[coll][z][e] = exRate1 + exRates2[coll][z][e]; // cm-3 s-1
            }
        }
    }

    // Copy exRates into exRateB, exRateC, exRateEF for convenience
    for (int z = 0; z < p.nbinsZ; ++z) {
        // float sumB = 0.0f;
        for (int e = 0; e < p.nbinsE; ++e) {
            exRateB[z][e]  = exRates[p.excitation_singlet_B][z][e];
            exRateC[z][e]  = exRates[p.excitation_singlet_C][z][e];
            exRateEF[z][e] = exRates[p.excitation_singlet_EF][z][e];
            // sumB += exRates[p.excitation_singlet_EF][z][e];
        }
        // std::cout << sumB << ", ";
    }
    // std::cout << "\n";

}

void Precip::sumQionOverE(){
        // Sum the ionisation rates over energy bins to get the total ionisation rate at each altitude
    for (int z = 0; z < p.nbinsZ; ++z) {
        for (int e = 0; e < p.nbinsE; ++e) {
            qz[z] += qion[e][z] * p.dE(e); // cm-1
            qz1[z] += qion1[e][z] * p.dE(e); // cm-1
        }
    }
}



void Precip::writeQionToFile(bool primariesOnly) {
    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "qione_" << suf << (primariesOnly ? "_primaries" : "") << ".dat";
    const std::string filename = ss.str();

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "qion[e][z]"},
        {"layout", "rows=e_index (0..nbinsE-1), cols=P_index (0..nbinsP-1)"},
        {"nbinsE", std::to_string(p.nbinsE)},
        {"Egrid_type", energyGrid.typestring()},
        {"Emin (eV)", std::to_string(p.E0)},
        {"Emax (eV)", std::to_string(p.Emax)},
        {"nbinsP", std::to_string(p.nbinsP)},
        {"Pgrid_type", "logarithmic"},
        {"P0 (Pa)", std::to_string(p.P0)},
        {"P1 (Pa)", std::to_string(p.P1)},
        {"units", "cm-1 eV-1"},
    };

    const bool ok = utils::io::write_dat_with_header(
        filename,
        "Precip model output: qion matrix",
        meta,
        /*column_names=*/{},  
        [&](std::ostream& os) {
            for (int e = 0; e < p.nbinsE; ++e) {
                for (int z = 0; z < p.nbinsZ; ++z) {
                    os << qion[e][z] << " ";
                }
                os << "\n";
            }
        }
    );

    if (ok) std::cout << "Wrote qion to " << filename << "\n";
    else    std::cerr << "Failed to write qion to " << filename << "\n";
}



void Precip::writeQzToFile(bool primariesOnly) {
    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "qion_" << suf << (primariesOnly ? "_primaries" : "") << ".dat";
    const std::string filename = ss.str();

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "qion"},
        {"units", " cm-1"},
        {"source_label", suf},
        {"nbinsP", std::to_string(p.nbinsP)},
        {"Pgrid_type", "logarithmic"},
        {"P0 (Pa)", std::to_string(p.P0)},
        {"P1 (Pa)", std::to_string(p.P1)},
    };

    const std::vector<std::string> cols = {"P [Pa]", "z [km]", "qion [cm-1]", "qion1 [cm-1]"};

    const bool ok = utils::io::write_dat_table_fixed_width(
        filename,
        "Precip model output: qion",
        meta,
        cols,
        p.nbinsZ,
        [&](int z, std::ostream& os, int w) {
            os << std::right << std::setw(w) << world.P[z]
            << " "        << std::setw(w) << world.Z[z]/1e3
            << " "        << std::setw(w) << qz[z]
            << " "        << std::setw(w) << qz1[z];
        },
        /*col_width=*/12,
        /*precision=*/6
    );

    if (ok) std::cout << "Wrote qz to " << filename << "\n";
}



void Precip::writeFUVExRatesToFile() {
    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "exratesFUV_" << suf << ".dat";
    const std::string filename = ss.str();

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "exRate[e][z]"},
        {"states", "B C EF"},
        {"layout",
         "rows=P_index (0..nbinsP-1). For each row, energies e=0..nbinsE-1 are written as triplets: "
         "[B(P,e) C(P,e) EF(P,e)] concatenated across e."},
        {"nbinsE", std::to_string(p.nbinsE)},
        {"Egrid_type", energyGrid.typestring()},
        {"Emin (eV)", std::to_string(p.E0)},
        {"Emax (eV)", std::to_string(p.Emax)},
        {"nbinsP", std::to_string(p.nbinsP)},
        {"Pgrid_type", "logarithmic"},
        {"P0 (Pa)", std::to_string(p.P0)},
        {"P1 (Pa)", std::to_string(p.P1)},
        {"units", "cm-1"},
    };

    const bool ok = utils::io::write_dat_with_header(
        filename,
        "Precip model output: excitation rates (B, C, EF)",
        meta,
        /*column_names=*/{}, // too many columns to name sensibly
        [&](std::ostream& os) {
            for (int z = 0; z < p.nbinsP; ++z) {
                for (int e = 0; e < p.nbinsE; ++e) {
                    os << exRates[p.excitation_singlet_B][z][e]  << " "
                       << exRates[p.excitation_singlet_C][z][e]  << " "
                       << exRates[p.excitation_singlet_EF][z][e] << " ";
                }
                os << "\n";
            }
        }
    );

    if (ok) std::cout << "Wrote FUV excitation rates to " << filename << "\n";
    else    std::cerr << "Failed to write FUV excitation rates to " << filename << "\n";
}

void Precip::writeExRatesToFile() {
    std::stringstream ss;
    const std::string suf = src->label(sp);
    ss << outdir << "exrates_" << suf << ".dat";
    const std::string filename = ss.str();

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "exRate[e][P]"},
        {"excitations", "all collisions"},
        {"layout",
         "rows=P_index (0..nbinsP-1). For each row, energies e=0..nbinsE-1 are written as lists: "
         "[elastic, ionisation, a, b, c, e, B, C, vibrational, rotational, EF] concatenated across e."},
        {"nbinsE", std::to_string(p.nbinsE)},
        {"Egrid_type", energyGrid.typestring()},
        {"Emin (eV)", std::to_string(p.E0)},
        {"Emax (eV)", std::to_string(p.Emax)},
        {"nbinsP", std::to_string(p.nbinsP)},
        {"Pgrid_type", "logarithmic"},
        {"P0 (Pa)", std::to_string(p.P0)},
        {"P1 (Pa)", std::to_string(p.P1)},
        {"units", "cm-1"},
    };

    const bool ok = utils::io::write_dat_with_header(
        filename,
        "Precip model output: excitation rates (all collisions)",
        meta,
        /*column_names=*/{}, // too many columns to name sensibly
        [&](std::ostream& os) {
            for (int z = 0; z < p.nbinsP; ++z) {
                for (int e = 0; e < p.nbinsE; ++e) {
                    for (int coll = 0; coll < p.nCollisions; ++coll) {
                        os << exRates[coll][z][e] << " ";
                    }
                }
                os << "\n";
            }
        }
    );

    if (ok) std::cout << "Wrote excitation rates to " << filename << "\n";
    else    std::cerr << "Failed to write excitation rates to " << filename << "\n";
}


void Precip::writePhiToFile(std::vector<std::vector<double>>& phi, const std::string& filename) {
    const std::string path = outdir + filename;

    std::vector<utils::io::MetaLine> meta = {
        {"run_ID", sp.runid},
        {"quantity", "phi[e][P]"},
        {"layout", "rows=e_index (0..nbinsE-1), cols=P_index (0..nbinsP-1)"},
        {"nbinsE", std::to_string(p.nbinsE)},
        {"Egrid_type", energyGrid.typestring()},
        {"Emin (eV)", std::to_string(p.E0)},
        {"Emax (eV)", std::to_string(p.Emax)},
        {"nbinsP", std::to_string(p.nbinsP)},
        {"Pgrid_type", "logarithmic"},
        {"P0 (Pa)", std::to_string(p.P0)},
        {"P1 (Pa)", std::to_string(p.P1)},
        {"units", "cm-2 s-1 sr-1 eV-1"}
    };

    const bool ok = utils::io::write_dat_with_header(
        path,
        "Precip model output: phi matrix",
        meta,
        /*column_names=*/{},
        [&](std::ostream& os) {
            const int nE = std::min<int>(p.nbinsE, static_cast<int>(phi.size()));
            for (int e = 0; e < nE; ++e) {
                const int nZ = std::min<int>(p.nbinsZ, static_cast<int>(phi[e].size()));
                for (int z = 0; z < nZ; ++z) {
                    os << phi[e][z] << " ";
                }
                os << "\n";
            }
        }
    );

    if (ok) std::cout << "Wrote phi to " << path << "\n";
    else    std::cerr << "Failed to write phi to " << path << "\n";
}


void Precip::run() {

    processPrimaries();
    nionToQion();
    processSecondaries(); 
    combinePrimarySecondary();
    writeQionToFile(false); // Write qion to a file, false means we write the full qion including secondaries
    sumQionOverE();
    writeQzToFile();
    writeFUVExRatesToFile();
    // writeExRatesToFile(); // Uncomment if you want to write all excitation rates to a file

}

