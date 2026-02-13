/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // IonDensity.hpp

#pragma once
#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include "World1D.hpp"
#include "Source.hpp"

// Ion density is computed by a chemistry model attached to the World.
// (Chosen per-world in WorldFactory.)


class IonDensity {
private:
    Params sp;                           // Simulation parameters
    std::shared_ptr<Source> src;         // Source of ionization data
    World1D &world;

    std::vector<double> Qtot;            // Total ionization rate at each altitude (cm^-3 s^-1)

    // Output directory (set in ctor using runid)
    std::string outdir;

    void computeIonDensity();
    void writeIonDensityToFile() const;

public:
    IonDensity(Params params_,
               std::shared_ptr<Source> src_,
               World1D &world_);

    void readIonisationRates(Params& params,
                             std::shared_ptr<Source> src,
                             const float F = 1.0f);

    void setIonisationRates(std::vector<float>& Qzin,
                            const float F = 1.0f);

    void run() {
        std::cout << "Computing ion density..." << std::endl;
        computeIonDensity();
        writeIonDensityToFile();
    }

    double getH3pColumnDensity() const;
};
