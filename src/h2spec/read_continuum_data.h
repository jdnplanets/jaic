/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines functions to read the continuum emission data from the MOLAT database:
 https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2cont97q5.php 
 Reference: Abgrall, H., Roueff, E., Liu, X., Shemansky D.E. 1997, Astroph. J., 481, 557-566
*/

 // read_continuum_data.h

#ifndef READ_CONTINUUM_DATA_H
#define READ_CONTINUUM_DATA_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <algorithm>

using namespace std;

struct ContinuumRecord {
    double ke;       // Energy value (from evfp)
    float a;         // First TRP column (from trp[ie][iif][0]) 
};

struct ContinuumData {
    int ev;  // The energy of the upper level
    vector<ContinuumRecord> records; // All emission records for this block (v value)
};

// Our map type: The outer key is the state label (ETAT), the inner key is a pair {vu, ju},
// and the value is a vector of ContinuumData blocks (one per block read from the file).
using ContinuumMap = map<string, map<pair<int,int>, map<int, ContinuumData>>>;

void loadContinuumFile(const std::string &filename, ContinuumMap &emap);
int getContinuumProbForLevel(const ContinuumMap &emap, const std::string &stateLabel, 
    int v, int ju, int jl, double& ev, vector<double>& ek, vector<double>& prob, bool verbose);
#endif // READ_CONTINUUM_DATA_H