/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */


 
 // Example usage of the JAIC code
 

 // main.cpp


#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>

#include "constants.h"
#include "WorldFactory.hpp"
#include "World1D.hpp"
#include "Sources.hpp"
#include "Precip.hpp"
#include "H2Column.hpp"
#include "H3pColumn.hpp"
#include "IonDensity.hpp"
#include "Conductivity.hpp"


using namespace std;


int ExampleJAICRun(); // Forward declaration of the MonoSingleRun function


int main() {

    ExampleJAICRun(); 

};


int ExampleJAICRun() {

    // Initialize the 1D world model
    auto world = make_world("Jupiter");



    float E = 30; // Energy in keV
    Params params{ 
        .N = 10000,             // Number of particles
        .Zinit = world->Z1,     // Initial position in m
        .Einit = E*1e3f,        // Initial energy in eV
        .runid = "example"      // Name   
    };

    // a mono‐energetic beam:
    auto src = make_shared<MonoSource>();
    // a BPL Isotropic beam:
    // E is ignored
    // auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
    //                                                       /*beta*/  10.0f, 
    //                                                       /*phi*/   15e3f);

    // a data source:
    // auto src = make_shared<DataSource>(/*filename*/"example_spectrum.txt");

    Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
    precip.run();

    double F = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1


    IonDensity iondens(params, src, *world);
    // iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    iondens.setIonisationRates(precip.qz, F); // ...or set the ionization rates from the values in precip
    iondens.run(); // Run the ion density calculation
    cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
    Conductivity conduct(params, src, *world);
    conduct.run(); // Run the conductivity calculation
    cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;

    H3pColumn h3pcol(params, src, *world);
    h3pcol.run();
    cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;

    H2Column h2col(params, src, *world);
    // h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
    h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F); // ...or set the excitation rates from the values in precip
    h2col.run(); // Run the H2 column calculations
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    return 0;
}

