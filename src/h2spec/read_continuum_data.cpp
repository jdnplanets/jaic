/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This is a fairly literal C++ port of the FORTRAN code to read the continuum emission data from the MOLAT database:
https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2cont97q5.php
Reference: Abgrall, H., Roueff, E., Liu, X., Shemansky D.E. 1997, Astroph. J., 481, 557-566
*/

 // read_continuum_data.cpp


#include "read_continuum_data.h"

using namespace std;

// Constants matching original Fortran dimensions.
const int IDIVE   = 100;    // For arrays EVE and SPO
const int IDIVF   = 2000;   // For arrays TRP, EVFP, EVFR
const int NUM_SPO = 5;      // Second dimension for SPO
const int NUM_TRP = 2;      // Third dimension for TRP

//---------------------------------------------------------------------
// Structures for storing dissociative–emission data.
//
// ContinuumRecord represents one “line” of dissociative emission output.
// (i.e. a record with IF-index, energy (ke), and two TRP values.)
//
// ContinuumData represents the data for one block. In our code we separate
// the emission data into two “blocks”: one for the first 50 records and one
// for the last 50 (or fewer) records, as in the original Fortran printing.
//---------------------------------------------------------------------


//---------------------------------------------------------------------
//  lcontrp:  Reads one “block” of data from the input stream.  
//           Each read operation first fetches an entire line via getline,
//           then uses an istringstream to parse the tokens.
//---------------------------------------------------------------------
void lcontrp(istream &ilec,
             int &ibe, int &iefin, int &iffip, int &iffir,
             vector<double> &evfp, vector<double> &evfr,
             string &etat,
             vector<double> &eve,
             vector<vector<float>> &spo,
             vector<vector<vector<float>>> &trp,
             int idive, int idivf)
{
    // Zero out arrays.
    for (int i = 0; i < idive; ++i) {
        eve[i] = 0.0;
        for (int j = 0; j < NUM_SPO; ++j) {
            spo[i][j] = 0.0f;
        }
        for (int j = 0; j < idivf; ++j) {
            trp[i][j][0] = 0.0f;
            trp[i][j][1] = 0.0f;
        }
    }

    string headerLine;
    do {
        getline(ilec, headerLine);
    } while (headerLine.empty());
    {
        istringstream headerStream(headerLine);
        headerStream >> ibe >> iefin >> iffip >> iffir >> etat;
    }

    // --------------------------
    // Read a dummy 80-character record
    // --------------------------
    {
        string dummy;
        getline(ilec, dummy);
    }

    int iemin = 1;
    int ifmin  = 1;
    
    // --------------------------
    // Read another dummy 80-character record.
    // --------------------------
    {
        string dummy;
        getline(ilec, dummy);
    }

    // --------------------------
    // Read the EVE array.
    // Fortran: READ(ILEC,370)(EVE(IE),IE=IEMIN,IEFIN) with FORMAT(1X,4(1X,F15.5))
    // We need to read (iefin - iemin + 1) double values. They may span several lines.
    // --------------------------
    int countEve = iefin - (iemin - 1);
    int numbersRead = 0;
    while (numbersRead < countEve) {
        string line;
        getline(ilec, line);
        if (line.empty())
            continue;
        istringstream iss(line);
        double value;
        while (iss >> value && numbersRead < countEve) {
            eve[iemin - 1 + numbersRead] = value;
            numbersRead++;
        }
    }

    // --------------------------
    // Read the SPO arrays for 4 states.
    // For each state, first a dummy 80-character record is read.
    // Then (iefin - iemin + 1) floats are read.
    // Fortran FORMAT: 37 = (1X,5E13.6)
    // --------------------------
    for (int state = 0; state < 4; ++state) {
        {
            string dummy;
            getline(ilec, dummy);
        }
        int countFloats = iefin - (iemin - 1);
        int floatsRead = 0;
        while (floatsRead < countFloats) {
            string line;
            getline(ilec, line);
            if (line.empty())
                continue;
            istringstream iss(line);
            float val;
            while (iss >> val && floatsRead < countFloats) {
                spo[iemin - 1 + floatsRead][state] = val;
                floatsRead++;
            }
        }
    }

    // --------------------------
    // Loop NDJ from 2 down to 1.
    // For each NDJ, set parameters and then loop IET (1 to 20) and IF.
    // For each IF, a line is read that contains one double and (iet2-iet1+1) floats.
    // Fortran FORMAT 304: (1X,F11.2,1X,4E13.6)
    // --------------------------
    for (int ndj = 2; ndj >= 1; --ndj) {
        int ibf, local_iefin;
        if (ndj == 2) {
            ibf = ibe - 1;
            local_iefin = iffir;
        } else { // ndj == 1
            ibf = ibe + 1;
            local_iefin = iffip;
        }
        
        {
            string dummy;
            getline(ilec, dummy);
        }
        
        for (int iet = 1; iet <= 20; ++iet) {
            int iet1 = iemin + 4 * (iet - 1);
            if (iet1 > iefin)
                break;  // As in Fortran: GO TO 223
            int iet2 = min(iemin + 4 * iet - 1, iefin);
            
            {
                string dummy;
                getline(ilec, dummy);
            }
            
            for (int iif = ifmin; iif <= local_iefin; ++iif) {
                string line;
                // Read the next nonempty line.
                do {
                    getline(ilec, line);
                } while (line.empty());

                istringstream iss(line);
                double temp;
                iss >> temp;
                if (ndj == 2)
                    evfr[iif - 1] = temp;
                else
                    evfp[iif - 1] = temp;
                // Now read the TRP values for IE from iet1 to iet2.
                for (int ie = iet1; ie <= iet2; ++ie) {
                    float val;
                    iss >> val;
                    trp[ie - 1][iif - 1][ndj - 1] = val;
                }
            }
                // --------------------------

        }
        // Read three final dummy 80-character records.
        // --------------------------
        for (int i = 0; i < 3; ++i) {
            string dummy;
            getline(ilec, dummy);
        }

    }

}


void lcontq(istream &ilec,
             int &ibe, int &ibf, int &iefin, int &iffiq,
             vector<double> &evfq,
             string &etat,
             vector<double> &eve,
             vector<vector<float>> &spo,
             vector<vector<float>> &trp,
             int idive, int idivf)
{
    // Zero out arrays.
    for (int i = 0; i < idive; ++i) {
        eve[i] = 0.0;
        for (int j = 0; j < NUM_SPO; ++j) {
            spo[i][j] = 0.0f;
        }
        for (int j = 0; j < idivf; ++j) {
            trp[i][j] = 0.0f;
        }
    }

    // --------------------------
    // Read the header line.
    // Fortran FORMAT: (1X,4I10,A4)
    // Expected tokens: IBE, IEFIN, IFFIP, IFFIR, ETAT
    // --------------------------
    string headerLine;
    do {
        getline(ilec, headerLine);
    } while (headerLine.empty());
    {
        istringstream headerStream(headerLine);
        headerStream >> ibe >> ibf >> iefin >> iffiq >> etat;
    }

    // --------------------------
    // Read a dummy 80-character record (as in Fortran READ(ILEC,63)...)
    // --------------------------
    {
        string dummy;
        getline(ilec, dummy);
    }

    int iemin = 1;
    int ifmin  = 1;
    
    // --------------------------
    // Read another dummy 80-character record.
    // --------------------------
    {
        string dummy;
        getline(ilec, dummy);
    }

    // --------------------------
    // Read the EVE array.
    // Fortran: READ(ILEC,370)(EVE(IE),IE=IEMIN,IEFIN) with FORMAT(1X,4(1X,F15.5))
    // We need to read (iefin - iemin + 1) double values. They may span several lines.
    // --------------------------
    int countEve = iefin - (iemin - 1);
    int numbersRead = 0;
    while (numbersRead < countEve) {
        string line;
        getline(ilec, line);
        if (line.empty())
            continue;
        istringstream iss(line);
        double value;
        while (iss >> value && numbersRead < countEve) {
            eve[iemin - 1 + numbersRead] = value;
            numbersRead++;
        }
    }

    // --------------------------
    // Read the SPO arrays for 4 states.
    // For each state, first a dummy 80-character record is read.
    // Then (iefin - iemin + 1) floats are read.
    // Fortran FORMAT: 37 = (1X,5E13.6)
    // --------------------------
    for (int state = 0; state < 2; ++state) {
        {
            string dummy;
            getline(ilec, dummy);
        }
        int countFloats = iefin - (iemin - 1);
        int floatsRead = 0;
        while (floatsRead < countFloats) {
            string line;
            getline(ilec, line);
            if (line.empty())
                continue;
            istringstream iss(line);
            float val;
            while (iss >> val && floatsRead < countFloats) {
                spo[iemin - 1 + floatsRead][state] = val;
                floatsRead++;
            }
        }
    }

    // --------------------------
    // Loop NDJ from 2 down to 1.
    // For each NDJ, set parameters and then loop IET (1 to 20) and IF.
    // For each IF, a line is read that contains one double and (iet2-iet1+1) floats.
    // Fortran FORMAT 304: (1X,F11.2,1X,4E13.6)
    // --------------------------

    int local_iefin;

    ibf = ibe - 1;
    local_iefin = iffiq;
    
    {
        string dummy;
        getline(ilec, dummy);
    }
    
    for (int iet = 1; iet <= 20; ++iet) {
        int iet1 = iemin + 4 * (iet - 1);
        if (iet1 > iefin)
            break;  // As in Fortran: GO TO 223
        int iet2 = min(iemin + 4 * iet - 1, iefin);
        
        {
            string dummy;
            getline(ilec, dummy);
        }
        
        for (int iif = ifmin; iif <= local_iefin; ++iif) {
            string line;
            // Read the next nonempty line.
            do {
                getline(ilec, line);
            } while (line.empty());

            istringstream iss(line);
            double temp;
            iss >> temp;

            evfq[iif - 1] = temp;
            // Now read the TRP values for IE from iet1 to iet2.
            for (int ie = iet1; ie <= iet2; ++ie) {
                float val;
                iss >> val;
                trp[ie - 1][iif - 1] = val;
            }
        }
            // --------------------------

    }
    // Read three final dummy 80-character records.
    // --------------------------
    for (int i = 0; i < 3; ++i) {
        string dummy;
        getline(ilec, dummy);
        // printf("Dummy5: %s\n", dummy.c_str());
    }
}


//---------------------------------------------------------------------
// loadContinuumFile: Opens a file, reads its data blocks (four per file),
// and appends each block's emission data into the provided map.
// Now we create one ContinuumData block per v value (v = ie-1).
// For each such block we record all IF records (from 1 to IFFIP).
// The map is keyed by state label (ETAT) and a pair {v, ju}.
//---------------------------------------------------------------------
void loadContinuumFile(const string &filename, ContinuumMap &emap)
{
    ifstream inFile("data/h2spec/abgrall_continuum/"+filename);
    if (!inFile) {
        cerr << "Error opening input file: " << filename << "\n";
        return;
    }

    // Scan through the file and count the number of occurrences of "E STATE HAS"
    int count = 0;
    string line;
    while (getline(inFile, line)) {
        size_t pos = 0;
        while ((pos = line.find("E S", pos)) != string::npos) {
            ++count;
            pos += string("E S").length();
        }
    }
    // Reset the file stream to the beginning for further processing
    inFile.clear();
    inFile.seekg(0, ios::beg);
    // cout << "Number of upper levels in file: " << count << "\n";
    
    // A simple test: if the filename contains "Q", assume Q branch.
    bool isQBranch = (filename.find("Q") != string::npos);
    bool isRPBranch = (filename.find("RP") != string::npos);
    bool isPBranch = false;
    if (!isQBranch && !isRPBranch) {
        isPBranch = true;
    }

    
    if (isQBranch || isPBranch) {
        // --- Process Q–branch file ---
        for (int block = 0; block < count; ++block) {
            int ibe, ibf, iefin, iffIQ;
            string etat;
            vector<double> EVE(IDIVE, 0.0);
            vector<double> evfq(IDIVF, 0.0);
            vector<vector<float>> SPO(IDIVE, vector<float>(5, 0.0f));
            vector<vector<float>> TRP(IDIVE, vector<float>(IDIVF, 0.0f));
            lcontq(inFile, ibe, ibf, iefin, iffIQ,
                   evfq, etat, EVE, SPO, TRP, IDIVE, IDIVF);
            // For Q branch, we take ju = IBF.
            int ju = ibf;
            int jl = ju;
            if (isPBranch) {
                ju = ibe - 1;
                jl = ju + 1;
            }
            int num_v = iefin;
            for (int v = 0; v < num_v; ++v) {
                ContinuumData ed;
                ed.ev = EVE[v];
                for (int iif = 0; iif < iffIQ; ++iif) {
                    ContinuumRecord rec;
                    rec.ke = evfq[iif];
                    rec.a = TRP[v][iif];  // Only one column in Q branch.
                    ed.records.push_back(rec);                }
                pair<int,int> keyPair = make_pair(v, ju);
                emap[etat][keyPair][jl] = ed;
            }
        }
    } else {
        // --- Process RP–branch file (as before) ---
        for (int block = 0; block < count; ++block) {
            vector<double> eve(IDIVE, 0.0);
            vector<double> evfp(IDIVF, 0.0);
            vector<double> evfr(IDIVF, 0.0);
            vector<vector<float>> spo(IDIVE, vector<float>(NUM_SPO, 0.0f));
            vector<vector<vector<float>>> trp(IDIVE, vector<vector<float>>(IDIVF, vector<float>(NUM_TRP, 0.0f)));
            int ibe, iefin, iffip, iffir;
            string etat;
            lcontrp(inFile, ibe, iefin, iffip, iffir,
                    evfp, evfr, etat, eve, spo, trp, IDIVE, IDIVF);
            int ju = ibe - 1;
            int jl = ju - 1;
            int num_v = iefin;
            for (int v = 0; v < num_v; ++v) {
                ContinuumData ed;
                ed.ev = eve[v];
                for (int iif = 0; iif < iffip; ++iif) {
                    ContinuumRecord rec;
                    rec.ke = evfp[iif];
                    rec.a = trp[v][iif][1];
                    ed.records.push_back(rec);
                }
                pair<int,int> keyPair = make_pair(v, ju);
                emap[etat][keyPair][jl] = ed;
            }
            jl = ju + 1;
            for (int v = 0; v < num_v; ++v) {
                ContinuumData ed;
                ed.ev = eve[v];
                for (int iif = 0; iif < iffip; ++iif) {
                    ContinuumRecord rec;
                    rec.ke = evfp[iif];
                    rec.a = trp[v][iif][0];
                    ed.records.push_back(rec);
                }
                pair<int,int> keyPair = make_pair(v, ju);
                emap[etat][keyPair][jl] = ed;
            }
        }
    }
}

//---------------------------------------------------------------------
int getContinuumProbForLevel(const ContinuumMap &emap,
                    const string &stateLabel,
                    int vu, int ju, int jl, double& ev, vector<double>& ek, vector<double>& prob, bool verbose)
{
    // Skip unphysical transitions.
    if (stateLabel == "B" && ju == jl) return 1;
    
    try {
        auto rec = emap.at(stateLabel).at({vu, ju}).at(jl);
        ev = rec.ev;
        for (const auto &r : rec.records) {
            ek.push_back(r.ke);
            prob.push_back(r.a);
        }
    } catch (const std::out_of_range& e) {
        if (verbose){
            cerr << "No continuum data for state=" << stateLabel
                    << ", v'=" << vu << ", J'=" << ju << ", J''='"
                    << jl << "\n";
        }
        return 1;
    }

    
    return 0;
}
