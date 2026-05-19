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
using namespace constants;

int MonoRunSingle(); // Forward declaration of the MonoSingleRun function
int MonoRunMultiConstNFlux(); // Forward declaration of the MonoMultiRun function
int MonoRunMultiConstEFlux(); // Forward declaration of the MonoMultiRun function with constant energy flux
int JunoRun(); // Forward declaration of the JunoRun function
int FACrun(); // Forward declaration of the FACRun function
int FACrunAur(); // Forward declaration of the FACRunAur function
int BPLGRunSingle(); // Forward declaration of the BPLRunSingle function
int BPLRunMulti(); // Forward declaration of the BPLRunMulti function
int h3p_h2_comparison_gridplot(); // Forward declaration of the h3p_h2_comparison_gridplot function
int FlareRunSingle(); // Forward declaration of the FlareRunSingle function

int FACrunAur_oneIndex(int idx); // Forward declaration of the FACrunAur_oneIndex function
int MonoRunMultiConstNFlux_oneIndex(int idx); // Forward declaration of the MonoMultiRun function with constant number flux for one index
int MonoRunMultiConstEFlux_oneIndex(int idx); // Forward declaration of the MonoMultiRun function with constant energy flux for one index
int BPLRunMulti_oneIndex(int idx); // Forward declaration of the BPLRunMulti function for one index

int main(int argc, char** argv) {

    cout << endl;
    cout << "*****************************************************************************" << endl;
    cout << "                  JAIC: Jupiter Auroral Ionosphere Code" << std::endl;
    cout << "*****************************************************************************" << endl;
    cout << endl;

    if (argc < 2) {
        MonoRunSingle(); 
        // MonoRunMultiConstNFlux();
        // MonoRunMultiConstEFlux();
        // JunoRun();
        // FACrun(); 
        // FACrunAur(); // Run the FAC run for aurora
        // BPLGRunSingle();
        // BPLRunMulti();
        // h3p_h2_comparison_gridplot();
        // FlareRunSingle();


    }
    else {
        int idx = std::atoi(argv[1]);
        // FACrunAur_oneIndex(idx);
        // MonoRunMultiConstNFlux_oneIndex(idx);
        MonoRunMultiConstEFlux_oneIndex(idx);
        // BPLRunMulti_oneIndex(idx);
    }
    return 0;
};


int MonoRunSingle() {

    // Initialize the 1D world model
    WorldOverrides o;
    o.Bdipang_deg = 90.0; // magnetic field dip angle in degrees (90 = vertical, 0 = horizontal)
    auto world = make_world("Jupiter", o);
    // world->write_species_densities("jupiter_species_densities_newp.dat"); // Write the species densities to a CSV file
    // exit(0);
    // WorldOverrides o;
    // o.g_m_s2 = 1000.;
    // o.B_T = 0.1;
    // auto world = make_world("Brown Dwarf", o);



    float E = 100; // Energy in keV
    Params params{ 
        .N = 10000,          // Number of particles
        .Zinit = world->Z1,   // Initial d position in m
        .Einit = E*1e3f,     // Initial energy in eV
        .runid = "colcounttest"      // Name   
    };

    // a mono‐energetic beam:
    auto src = make_shared<MonoSource>();
    Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
    precip.run();
    // precip.runPrimariesOnly(); // Process the primary electrons only
    // precip.runSecondariesOnly();

    double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
    // double Ef = 100; // Energy flux in mW m-2
    // double F = Ef * 1e-3*1e-4/(E * 1e3 * qe); // Number flux to carry Ef mW m-2 flux in cm-2 s-1
    double F = F0; // Number flux to carry E mW m-2 flux in cm-2 s-1


    IonDensity iondens(params, src, *world);
    // iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
    iondens.run(); // Run the ion density calculation
    cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
    Conductivity conduct(params, src, *world);
    conduct.run(); // Run the conductivity calculation
    cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;


    // H2Column h2col(params, src, *world);
    // h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
    // // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
    // h2col.run(); // Run the H2 column calculations
    // cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    // cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    // cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    // H3pColumn h3pcol(params, src, *world, true); // true to apply non-LTE scaling factors
    // h3pcol.run();
    // cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;
    // int peakzIndex = h3pcol.getPeakzIndex();
    // cout << "Peak H3+ emission altitude: " << world->Z[peakzIndex]/1e3 << " km, where T = " << world->T[peakzIndex] << " K" << endl;


    return 0;
}


int MonoRunMultiConstNFlux() {


    auto world = make_world("Jupiter"); // Initialize the 1D world model

        // // Create a log-spaced array of 100 points between 10 and 1000

    vector<float> Es;
    int n_points = 100;
    float start = 1.0f;
    float end = 1000.0f;
    Es.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
        Es.push_back(val);
    }
    
    int j = -1;
    for (float E : Es) {
        j++;
        if (j < 91) continue;
        cout << "Running for E[" << j << "] = " << E << " keV" << endl;

        Params params{ 
            .N = 100000,          // Number of particles
            .Zinit = world->Z1,       // Initial position in m
            .Einit = E*1e3f,        // Initial energy in eV
            .runid = "MonoE_const_nFlux_1uAm2_v5" // Name for the run
        };

        // a mono‐energetic beam:
        auto src = make_shared<MonoSource>();
        Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
        precip.runlite();

        //constant number flux
        double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
        double F = F0;
        // double Ef = 1; // Energy flux in mW m-2
        // double F = Ef * 1e-3*1e-4/(E * 1e3 * qe); // Number flux to carry Ef mW m-2 flux in cm-2 s-1. For Benmahi plot


        IonDensity iondens(params, src, *world); // 10 keV
        iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
        // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
        iondens.run(); // Run the ion density calculation
        Conductivity conduct(params, src, *world);
        conduct.run(); // Run the conductivity calculation


        // H2Column h2col(params, src, *world);
        // h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
        // // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
        // h2col.run(); // Run the H2 column calculations
        // cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
        // cout << "Total intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
        // cout << "Colour ratio: " << h2col.getColourRatio() << endl;

        H3pColumn h3pcol(params, src, *world, true);
        h3pcol.run();
        cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;

    }

    return 0;
}


// This is for the Benmahi Fig 8 comparison (with Ef=1 mW m-2) plot and SigmaP versus E plots
int MonoRunMultiConstEFlux() {

    auto world = make_world("Jupiter"); // Initialize the 1D world model

        // // Create a log-spaced array of 100 points between 10 and 1000

    float Efs[] = {1, 0.1, 10, 100}; // All for SigmaP versus E plot
    // float Efs[] = {1};
    for (float Ef : Efs) {

        vector<float> Es;
        int n_points = 100;
        float start = 1.0f; // keV
        float end = 1000.0f;
        Es.reserve(n_points);
        for (int i = 0; i < n_points; ++i) {
            float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
            Es.push_back(val);
        }
        
        int j = -1;
        for (float E : Es) {
            j++;
            // if (j < 96) continue; // Useful to restart if the GPU session is timed out
            cout << "Running for E[" << j << "] = " << E << " keV" << endl;
            string runid;
            if (Ef < 1) {
                std::ostringstream oss;
                oss << "MonoE_constEFlux_" << std::fixed << std::setprecision(1) << Ef << "mWm2_v5";
                runid = oss.str();
            }
            else {
                runid = "MonoE_constEFlux_" + std::to_string(static_cast<int>(Ef)) + "mWm2_v5";
            }

            Params params{ 
                .N = 100000,          // Number of particles
                .Zinit = world->Z1,       // Initial position in m
                .Einit = E*1e3f,        // Initial energy in eV
                .runid = runid // Name for the run
            };

            // a mono‐energetic beam:
            auto src = make_shared<MonoSource>();
            // Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
            // precip.runlite();

            // //constant number flux
            // // double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
            // // double F = F0;
            // double Ef = 1; // Energy flux in mW m-2
            double F = Ef * 1e-3*1e-4/(E * 1e3 * qe); // Number flux to carry Ef mW m-2 flux in cm-2 s-1. For Benmahi plot
            // cout << F << endl;
            // exit(0);

            IonDensity iondens(params, src, *world); // 10 keV
            // params.runid = "MonoE_constEFlux_1mWm2";
            iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
            // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
            iondens.run(); // Run the ion density calculation
            params.runid = runid; // Restore the original runid
            Conductivity conduct(params, src, *world);
            conduct.run(); // Run the conductivity calculation
            // exit(0);

            // H2Column h2col(params, src, *world);
            // // params.runid = "MonoE_constEFlux_1mWm2";
            // h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
            // // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
            // h2col.run(); // Run the H2 column calculations
            // cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
            // cout << "Total intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
            // cout << "Colour ratio: " << h2col.getColourRatio() << endl;
            
            H3pColumn h3pcol(params, src, *world, true);
            h3pcol.run();
            cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;

            

        }
    }

    return 0;
}



int BPLGRunSingle(){


    auto world = make_world("Jupiter"); // Initialize the 1D world model

    double jpi = 0.2; // field-aligned current density in μA m^-2
    double wth = 2.5e3; // Thermal energy in eV
    double n = 0.02 * 1e6; // Density in m^-3
    double jpi0 = e * n * sqrt(wth * eV / (2 * pi * m_e));
    float phi = wth * jpi*1e-6 / jpi0; // Knight voltage in eV for the given jpi and thermal energy, using the Knight relation

    float E = phi/1e3; // Energy in keV
    Params params{ 
        .N = 100000,          // Number of particles
        .Zinit = world->Z1,   // Initial position in m
        .Einit = E*1e3f,     // Initial energy in eV
        .runid = "stan_facs"      // Name   
    };

    // a BPL Iso beam:
    // auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
    //                                                     /*beta*/  10.0f, 
    //                                                     /*phi*/   phi);
    auto src = make_shared<BrokenPowerLawGaussianIsotropicSource>(/*alpha*/ 2.0f, 
                                                            /*beta*/  8.0f, 
                                                            /*phi*/   phi,
                                                            /*sigma*/ phi * (pow(10.0, 0.05) - 1),
                                                            /*peakFactor*/ 8.0f);
    Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
    precip.run();

    double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
    double F = F0 * jpi; // Number flux to carry jpi μA m-2 flux in cm-2 s-1

    IonDensity iondens(params, src, *world);
    iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
    iondens.run(); // Run the ion density calculation
    cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
    Conductivity conduct(params, src, *world);
    conduct.run(); // Run the conductivity calculation
    cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;

    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
    // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
    h2col.run(); // Run the H2 column calculations
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    H3pColumn h3pcol(params, src, *world, true);
    h3pcol.run();
    cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;


    return 0;

}


int JunoRun() {

    auto world = make_world("Jupiter"); // Initialize the 1D world model

    Params params{ 
        .N = 100000,          // Number of particles
        .Zinit = world->Z1,       // Initial position in m
        .Einit = 1,        // Initial energy in eV. Not used with this source
        .runid = "Ebert2019_4b" // Name for the run
        // .runid = "Ebert2019_4d" // Name for the run
        // .runid = "Ebert2021_6c" // Name for the run

    };

    
    auto src = make_shared<DataSource>(/*filename*/params.runid + ".txt");

    Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
    // precip.runlite();

    double F = src->getNumberFlux(); // Number flux in cm-2 s-1 from DataSource
    cout << "F = " << F << " cm-2 s-1" << endl;

    IonDensity iondens(params, src, *world);
    iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
    iondens.run(); // Run the ion density calculation
    cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
    Conductivity conduct(params, src, *world);
    conduct.run(); // Run the conductivity calculation
    cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;


    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
    // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
    h2col.run(); // Run the H2 column calculations
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    H3pColumn h3pcol(params, src, *world, true);
    h3pcol.run();
    cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;


    return 0;
}

// For djpi = 0.01 to 1.5 in steps of 0.01
// CPU Time:	10:01:26
// Maximum memory usage per node:	459.16MiB
// This run is for the  SigmaP versus Jpari results
int FACrun() {

    auto world = make_world("Jupiter"); // Initialize the 1D world model

    double djs[] = {0.5, 0.01};
    int pass = -1;
    for (double djpi : djs) {
        pass++;

        float ns[] = {0.02e6, 0.001e6, 0.05e6}; // All for SigmaP versus E plot
        // float ns[] = {0.001e6}; // All for SigmaP versus E plot

        for (float n : ns) {

            cout << "Processing for n = " << n/1e6 << " cm^-3" << endl;

            // double djpi = 0.01; // Step size for jpi
            cout << "Processing energies for jpi = " << djpi << " to 1.5 in steps of " << djpi << endl;
            double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1

            double j0, j1;
            if (n == 0.001e6) {
                j0 = 0.01;
                j1 = 0.5;
            }
            else if (n == 0.02e6) {
                j0 = 0.01;
                j1 = 3.0;
            }
            else if (n == 0.05e6) {
                j0 = 0.01;
                j1 = 3.0;
            }
            // j0 = 0.01;
            // j1 = 3.0;
            for (double jpi = j0; jpi <= j1; jpi += djpi) {
            // for (double jpi = 0.08; jpi <= 0.5; jpi += djpi) {

                double wth = 2.5e3; // Thermal energy in keV
                // double n = 0.015e6; // Density in cm^-3
                double jpi0 = e * n * sqrt(wth * eV / (2 * pi * m_e));
                float phi = wth * jpi*1e-6 / jpi0;

                cout << jpi0 << " uAm-2: ";
                cout << "j = " << jpi << " uAm-2, E = " << phi/1e3 << " keV" << endl;

                int nlab = n / 1e3; // For runid
                char nlab_buf[4];
                std::snprintf(nlab_buf, sizeof(nlab_buf), "%03d", nlab);
                std::string runid;
                if (pass == 0) {
                    runid = std::string("SigmaP_v_jpari_BPLG_n") + nlab_buf + "_coarse_v5";
                    cout << "Run ID: " << runid << endl;
                }
                else {
                    runid = std::string("SigmaP_v_jpari_BPLG_n") + nlab_buf+ "_v5";
                    cout << "Run ID: " << runid << endl;
                }
                // exit(0);
            
                Params params{ 
                    .N = 100000,                  // Number of particles
                    .Zinit = 3000e3f,             // Initial position in m
                    .Einit = phi,                 // Initial energy in eV
                    .runid = runid     // Name
                };

                // auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
                //                                                     /*beta*/  10.0f, 
                //                                                     /*phi*/   phi);
                auto src = make_shared<BrokenPowerLawGaussianIsotropicSource>(/*alpha*/ 2.0f, 
                                                                        /*beta*/  8.0f, 
                                                                        /*phi*/   phi,
                                                                        /*sigma*/ phi * (pow(10.0, 0.05) - 1),
                                                                        /*peakFactor*/ 8.0f);
 
                // auto src = make_shared<GaussianSource>(/*sigma*/ 0.1f*phi);

                Precip precip(params, src, *world);
                precip.runlite();
                // // exit(0);

                double F = F0 * jpi; // Number flux to carry the desired current in cm-2 s-1

                IonDensity iondens(params, src, *world); // 10 keV
                iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
                // iondens.setIonisationRates(precip.qz, F); // Set the
                iondens.run(); // Run the ion density calculation
                Conductivity conduct(params, src, *world);
                conduct.run(); // Run the conductivity calculation

                H3pColumn h3pcol(params, src, *world, true);
                h3pcol.run();


            }
            cout << "Finished processing all energies." << endl;

        }
        cout << "Finished processing all densities." << endl;

    }
    cout << "Finished processing all current steps" << endl;
   
    return 0;
}

// This run is for the intensity versus CR results
int FACrunAur() {

    auto world = make_world("Jupiter"); // Initialize the 1D world model

    double djpi = 0.01; // Step size for jpi
    cout << "Processing energies for jpi = 0.1 to 5 in steps of " << djpi << endl;
    double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1

    int j = -1;
    for (double jpi = djpi; jpi <= 5; jpi += djpi) {
    // for (double jpi = 4.5; jpi <= 5; jpi += djpi) {
        j++;
        // if (j < 32) continue;

        double wth = 2.5e3; // Thermal energy in keV
        double n = 0.02e6; // Density in cm^-3
        double jpi0 = e * n * sqrt(wth * eV / (2 * pi * m_e));
        float phi = wth * jpi*1e-6 / jpi0;

        cout << jpi0 << " uAm-2: ";
        cout << "j = " << j << " jpari = " << jpi << " uAm-2, E = " << phi/1e3 << " keV" << endl;
    
        int nlab = n / 1e3; // For runid
        char nlab_buf[4];
        std::snprintf(nlab_buf, sizeof(nlab_buf), "%03d", nlab);
        std::string runid = std::string("SigmaP_v_jpari_BPLG_n") + nlab_buf;
        cout << "Run ID: " << runid << endl;

        Params params{ 
            .N = 100000,                  // Number of particles
            .Zinit = 3000e3f,             // Initial position in m
            .Einit = phi,                 // Initial energy in eV
            .runid = runid
        };

        // // a mono‐energetic beam:
        // auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
        //                                                       /*beta*/  10.0f, 
        //                                                       /*phi*/   phi);


        auto src = make_shared<BrokenPowerLawGaussianIsotropicSource>(/*alpha*/ 2.0f, 
                                                            /*beta*/  8.0f, 
                                                            /*phi*/   phi,
                                                            /*sigma*/ phi * (pow(10.0, 0.05) - 1),
                                                            /*peakFactor*/ 8.0f);
        double F = F0 * jpi; // Number flux to carry the desired current in cm-2 s-1


        H2Column h2col(params, src, *world);
        // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
        h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
        h2col.run(); // Run the H2 column calculations
        cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
        cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
        cout << "Colour ratio: " << h2col.getColourRatio() << endl;

 
 
        IonDensity iondens(params, src, *world); // 10 keV
        iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
        iondens.run(); // Run the ion density calculation

        H3pColumn h3pcol(params, src, *world, true);
        h3pcol.run();
        cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;


        // #include <fstream>
        // std::ofstream outfile("maxj.txt");
        // outfile << j << std::endl;
        // outfile.close();



    }
    cout << "Finished processing all currents." << endl;
   
    return 0;
}





int BPLRunMulti() {


    auto world = make_world("Jupiter"); // Initialize the 1D world model

        // // Create a log-spaced array of 100 points between 10 and 1000

    vector<float> Es;
    int n_points = 100;
    float start = 0.1f;
    float end = 1000.0f;
    Es.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
        Es.push_back(val);
    }

    // vector<float> Es = {1, 10, 100, 1000};
    // vector<float> Es = {1000};

    
    int j = -1;
    for (float E : Es) {
        j++;
        // if (j < 96) continue;
        cout << "Running for E[" << j << "] = " << E << " keV" << endl;

        Params params{ 
            .N = 100000,          // Number of particles
            .Zinit = world->Z1,       // Initial position in m
            .Einit = E*1e3f,        // Initial energy in eV
            .runid = "BPL_const_nFlux_1uAm2" // Name for the run
        };

        // a mono‐energetic beam:
        auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
                                                              /*beta*/  10.0f, 
                                                              /*phi*/   E*1e3f);

        Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
        precip.runlite();


        double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
        double F = 1 * F0;

        IonDensity iondens(params, src, *world);
        iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
        // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
        iondens.run(); // Run the ion density calculation
        cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
        Conductivity conduct(params, src, *world);
        conduct.run(); // Run the conductivity calculation
        cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;


        // H2Column h2col(params, src, *world);
        // h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
        // // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
        // h2col.run(); // Run the H2 column calculations
        // cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
        // cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
        // cout << "Colour ratio: " << h2col.getColourRatio() << endl;

        H3pColumn h3pcol(params, src, *world, true);
        h3pcol.run();
        cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;

    }

    return 0;
}


// This is for the H3+ radiance comparison plot with H2 column results
int h3p_h2_comparison_gridplot() {


    auto world = make_world("Jupiter"); // Initialize the 1D world model

     // Create a log-spaced array of 100 points between 0.01 and 10
    vector<float> js;
    int n_points_j = 100;
    float start = 0.01f; // uA m-2
    float end = 10.0f;   // uA m-2
    js.reserve(n_points_j);
    for (int i = 0; i < n_points_j ; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points_j - 1));
        js.push_back(val);
    }
    
    // Create a log-spaced array of 100 points between 1 and 1000
    vector<float> Es;
    int n_points_E = 100;
    start = 1.0f;
    end = 1000.0f;
    Es.reserve(n_points_E);
    for (int i = 0; i < n_points_E; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points_E - 1));
        Es.push_back(val);
    }

    vector<vector<double>> H3p_radiances = vector<vector<double>>(n_points_j, vector<double>(n_points_E, 0.0));
    
    int j = -1;
    for (float jpi : js) {
        j++;
        int e = -1;
        for (float E : Es) {
            e++;

            cout << "Running for  jpi[" << j << "] = " << jpi << " uA m^-2, E[" << e << "] = " << E << " keV" << endl;
            Params params{ 
                .N = 100000,          // Number of particles
                .Zinit = world->Z1,       // Initial position in m
                .Einit = E*1e3f,        // Initial energy in eV
                .runid = "MonoE_const_nFlux_1uAm2_v5" // Name for the run
            };

            // a mono‐energetic beam:
            auto src = make_shared<MonoSource>();
            // auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
            //                                                       /*beta*/  10.0f, 
            //                                                       /*phi*/   E*1e3f);

            //constant number flux
            double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
            double F = jpi * F0;


            IonDensity iondens(params, src, *world); // 10 keV
            iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
            // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
            // Rename runid for this comparison plot
            params.runid = "H3P_H2_comparison_plot_v5";
            iondens.run(); // Run the ion density calculation
            Conductivity conduct(params, src, *world);
            conduct.run(); // Run the conductivity calculation

            H3pColumn h3pcol(params, src, *world, true);
            h3pcol.run();
            H3p_radiances[j][e] = h3pcol.getF335MRadiance();
            cout << "F335M Radiance: " << H3p_radiances[j][e] / 1e-6 << " μW m^-2 sr^-1" << endl;
            // exit(0);

        }

    }

    string filename = "out/H3p_radiances_nonLTE_v5.txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return 1;
    }

    // Write data
    for (int j = 0; j < n_points_j; j++) {
        for (int e = 0; e < n_points_E; e++) {
            outfile << "\t" << H3p_radiances[j][e];
        }
        outfile << "\n";
    }
    outfile.close();
    cout << "H3p radiances written to " << filename << endl;
    

    return 0;
}




int FlareRunSingle() {

    // Initialize the 1D world model
    auto world = make_world("Jupiter");
    


    float E = 30; // Energy in keV
    Params params{ 
        .N = 1000,          // Number of particles
        .Zinit = world->Z1,   // Initial position in m
        .Einit = E*1e3f,     // Initial energy in eV
        .runid = "flare"      // Name   
    };

    // a mono‐energetic beam:
    auto src = make_shared<MonoSource>();
    Precip precip(params, src, *world); // Initialize the precipitation simulation with the world model and source
    precip.run();

    double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
    // double Ef = 100; // Energy flux in mW m-2
    // double F = Ef * 1e-3*1e-4/(E * 1e3 * qe); // Number flux to carry Ef mW m-2 flux in cm-2 s-1
    double F = F0; // Number flux to carry E mW m-2 flux in cm-2 s-1


    // IonDensity iondens(params, src, *world);
    // // iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    // iondens.setIonisationRates(precip.qz, F); // Set the ionization rates based on the mixing ratios and the total ionization rate
    // iondens.run(); // Run the ion density calculation
    // cout << "H3+ column density: " << iondens.getH3pColumnDensity() << " cm^-2" << endl;
    // Conductivity conduct(params, src, *world);
    // conduct.run(); // Run the conductivity calculation
    // cout << "Pedersen conductance: " << conduct.getConductance() << " mho" << endl;


    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
    // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
    h2col.run(); // Run the H2 column calculations
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    // H3pColumn h3pcol(params, src, *world, true);
    // h3pcol.run();
    // cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;
    // int peakzIndex = h3pcol.getPeakzIndex();
    // cout << "Peak H3+ emission altitude: " << world->Z[peakzIndex]/1e3 << " km, where T = " << world->T[peakzIndex] << " K" << endl;


    return 0;
}


//########################################################################
// Setup for array job


int FACrunAur_oneIndex(int idx) {

    // Scan definition
    const double djpi   = 0.01;
    const double jpiMin = 1.0;   // 0.01
    const double jpiMax = 2.0;   // 3.0

    const int N = (int)std::llround((jpiMax - jpiMin) / djpi) + 1; // 500

    if (idx < 1 || idx > N) {
        cerr << "Error: idx=" << idx << " out of range 1.." << N << "\n";
        return 2;
    }

    const int j = idx - 1;
    const double jpi = jpiMin + j * djpi;

    cout << "Running idx=" << idx << " (j=" << j << "), jpi=" << jpi << " uA m-2\n";

    auto world = make_world("Jupiter"); // Initialize the 1D world model

    const double F0 = 6.25e8; // cm-2 s-1 for 1 uA m-2

    double wth = 2.5e3;   // keV
    double n   = 0.02e6;  // cm^-3

    double jpi0 = e * n * sqrt(wth * eV / (2 * pi * m_e));
    float phi  = (float)(wth * jpi * 1e-6 / jpi0); // eV

    cout << jpi0 << " uAm-2: "
         << "j = " << j << " jpari = " << jpi << " uAm-2, E = " << phi/1e3 << " keV\n";

    // Constant runid (shared directory for all tasks)
    int nlab = (int)(n / 1e3);
    char nlab_buf[4];
    std::snprintf(nlab_buf, sizeof(nlab_buf), "%03d", nlab);
    std::string runid = std::string("SigmaP_v_jpari_BPLG_n") + nlab_buf;

    cout << "Run ID: " << runid << "\n";

    Params params{
        .N     = 100000,
        .Zinit = 3000e3f,
        .Einit = phi,
        .runid = runid
    };

    auto src = make_shared<BrokenPowerLawGaussianIsotropicSource>(/*alpha*/ 2.0f, 
                                                            /*beta*/  8.0f, 
                                                            /*phi*/   phi,
                                                            /*sigma*/ phi * (pow(10.0, 0.05) - 1),
                                                            /*peakFactor*/ 8.0f);
    double F = F0 * jpi;

    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F);
    h2col.run();
    // cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR\n";
    // cout << "Total observed intensity: "   << h2col.getTotalObservedIntensity()   << " kR\n";
    // cout << "Colour ratio: "               << h2col.getColourRatio()              << "\n";

    // IonDensity iondens(params, src, *world); // 10 keV
    // iondens.readIonisationRates(params, src, F); // Read the ionization rates from a file
    // iondens.run(); // Run the ion density calculation

    // H3pColumn h3pcol(params, src, *world, true);
    // h3pcol.run();
    // cout << "Total H3+ radiance: " << h3pcol.getTotalRadiance() / 1e-6 << " μW m^-2 sr^-1" << endl;


    return 0;
}



int MonoRunMultiConstNFlux_oneIndex(int idx) {


    auto world = make_world("Jupiter"); // Initialize the 1D world model

        // // Create a log-spaced array of 100 points between 10 and 1000

    vector<float> Es;
    int n_points = 100;
    float start = 1.0f;
    float end = 1000.0f;
    Es.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
        Es.push_back(val);
    }

    // validate idx
    if (idx < 1 || idx > n_points) {
        cerr << "Error: idx=" << idx << " out of range 1.." << n_points << "\n";
        return 2;
    }

    int j = idx - 1;
    float E = Es[j];
    cout << "Running for E[" << j << "] = " << E << " keV" << endl;

    Params params{
        .N = 100000,                      // Number of particles
        .Zinit = world->Z1,               // Initial position in m
        .Einit = E * 1e3f,                // Initial energy in eV
        .runid = "MonoE_const_nFlux_1uAm2_v5"
    };

    auto src = make_shared<MonoSource>();

    // constant number flux
    double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
    double F = F0;


    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F);
    // params.runid += "_grodenT";
    h2col.run();
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    return 0;
}


int MonoRunMultiConstEFlux_oneIndex(int idx) {

    auto world = make_world("Jupiter"); // Initialize the 1D world model

    const double Ef = 1.0; // mW m-2

    vector<float> Es;
    int n_points = 100;
    float start = 1.0f; // keV
    float end = 1000.0f;
    Es.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
        Es.push_back(val);
    }

    // validate idx
    if (idx < 1 || idx > n_points) {
        cerr << "Error: idx=" << idx << " out of range 1.." << n_points << "\n";
        return 2;
    }

    int j = idx - 1;
    float E = Es[j];
    cout << "Running for E[" << j << "] = " << E << " keV" << endl;

    char runid_buf[64];
    if (Ef < 1.0) {
        std::snprintf(runid_buf, sizeof(runid_buf), "MonoE_constEFlux_%.1fmWm2_v5", Ef);
    } else {
        std::snprintf(runid_buf, sizeof(runid_buf), "MonoE_constEFlux_%dmWm2_v5", static_cast<int>(Ef));
    }
    std::string runid(runid_buf);

    Params params{
        .N = 100000,                      // Number of particles
        .Zinit = world->Z1,               // Initial position in m
        .Einit = E * 1e3f,                // Initial energy in eV
        .runid = runid
    };

    auto src = make_shared<MonoSource>();

    double F = Ef * 1e-3 * 1e-4 / (E * 1e3 * qe); // Number flux to carry Ef mW m-2 flux in cm-2 s-1

    H2Column h2col(params, src, *world);
    h2col.readExcitationRates(params, src, F);
    h2col.run();
    cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
    cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
    cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    return 0;
}


int BPLRunMulti_oneIndex(int idx) {


    auto world = make_world("Jupiter"); // Initialize the 1D world model

        // // Create a log-spaced array of 100 points between 10 and 1000

    vector<float> Es;
    int n_points = 100;
    float start = 1.0f;
    float end = 1000.0f;
    Es.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        float val = start * std::pow(end / start, static_cast<float>(i) / (n_points - 1));
        Es.push_back(val);
    }   

    // validate idx
    if (idx < 1 || idx > n_points) {
        cerr << "Error: idx=" << idx << " out of range 1.." << n_points << "\n";
        return 2;
    }   

    int j = idx - 1;
    float E = Es[j];
    cout << "Running for E[" << j << "] = " << E << " keV" << endl;
    
        Params params{ 
            .N = 100000,          // Number of particles
            .Zinit = world->Z1,       // Initial position in m
            .Einit = E*1e3f,        // Initial energy in eV
            .runid = "BPL_const_nFlux_1uAm2" // Name for the run
        };

        // a mono‐energetic beam:
        auto src = make_shared<BrokenPowerLawIsotropicSource>(/*alpha*/ 3.0f, 
                                                              /*beta*/  10.0f, 
                                                              /*phi*/   E*1e3f);

            // constant number flux
        double F0 = 6.25e8; // Number flux for 1 uA m-2, in cm-2 s-1
        double F = F0;
        H2Column h2col(params, src, *world);
        h2col.readExcitationRates(params, src, F); // Read the excitation rates from a file
        // h2col.setExcitationRates(precip.exRateB, precip.exRateC, precip.exRateEF, F);
        h2col.run(); // Run the H2 column calculations
        cout << "Total unabsorbed intensity: " << h2col.getTotalUnabsorbedIntensity() << " kR" << endl;
        cout << "Total observed intensity: " << h2col.getTotalObservedIntensity() << " kR" << endl;
        cout << "Colour ratio: " << h2col.getColourRatio() << endl;

    return 0;
}