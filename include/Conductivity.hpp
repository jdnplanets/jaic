/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // This class computes the conductivity based on ionization rates and species properties.

 // Conductivity.hpp

#pragma once
#include <vector>
#include <memory>
#include "World1D.hpp"
#include "SpeciesData.hpp"
#include "Source.hpp"


class Conductivity {

    std::vector<double> sigmaP_H3p;     // Conductivity values in mho
    std::vector<double> sigmaP_CH5p;    
    std::vector<double> sigmaP_C3Hnp;
    std::vector<double> sigmaP_e;       // Electron conductivity
    std::vector<double> sigmaP;         // Total Pedersen conductivity (sum of all species)


    static inline const std::string jaicroot = std::string(getenv("JAIC_ROOT")) + "/";
    std::string outdir = jaicroot + "out/conduct/";

private:
    World1D &world;
    Params sp;       
    std::shared_ptr<Source> src;

    
    void computeConductivity();
    void writeConductivityToFile() const;


public:
    // Constructor receives world with ion densities
    Conductivity(Params params_, 
                 std::shared_ptr<Source> src_, 
                 World1D &world_);

    void run() {
        std::cout << "Computing Pedersen conductivities..." << std::endl;
        computeConductivity();      // Compute conductivity based on ion density
        writeConductivityToFile(); // Write conductivity to file
    }

    double getConductance();

};
