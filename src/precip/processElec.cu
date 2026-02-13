/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file contains the main processing loop for the primary electrons. It advances each electron, 
checks for collisions, and processes any collisions that occur. It also keeps track of how many 
electrons are still active (i.e. not lost or dead) and how many are lost through each loss mechanism 
(top, bottom, sides, energy). The main kernel is processElectronsGPU(), which is called from 
runPrimariesKernel() in main.cpp. 

NOTE: The code contains conditions on __CUDA_ARCH__ to allow the file to be compiled if the project
is to be built without GPU support (e.g. if only later modules need to be run). The code DOES NOT RUN
if compiled without GPU support, but it compiles successfully and allows other modules to be built and run. 
*/

 // processElec.cu

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "Rnd.h"
#include "constants.h"
#include "simparams.h"
#include "device_arrays.h"
#include "Egrid.hpp"
#include "cuda_runtime.h"
#include <curand_kernel.h>
using namespace std;
using namespace constants;

#define CUDA_CHECK(call)                                                   \
    do {                                                                    \
        cudaError_t err = call;                                             \
        if (err != cudaSuccess) {                                           \
            fprintf(stderr,                                                 \
                    "CUDA error in %s (%s:%d): %s\n",                       \
                    #call, __FILE__, __LINE__, cudaGetErrorString(err));   \
            exit(err);                                                      \
        }                                                                   \
    } while (0)



__device__ curandState* gStates;

// getRandom:  use curand_uniform() when compiled for device;  rnd() when for host.
__host__ __device__
inline float getRandom(int tid)
{
#ifdef __CUDA_ARCH__
    // Device side → use cuRAND.  Caller must supply a valid tid.
    return curand_uniform(&gStates[tid]);
#else
    // Host side → use the CPU RNG.
    (void) tid;
    return rnd();
#endif
}

struct ElectronLossCounters {
    int Nactive;
    int Nlost_top;
    int Nlost_bottom;
    int Nlost_sides;
    int Nlost_energy;
};

__device__ float Elost_primaries = 0.0f, 
                 Edied_primaries = 0.0f, 
                 Egained_secondaries = 0.0f,
                 Ealive = 0.0f;


__device__ __host__
bool test_collision(const float ndens, 
                    const float total_sigma,
                    DeviceArrays a, 
                    const int i) {

    // float v = sqrt(2 * E * eV / m_e); // m/s
    float gamma = 1.0 + a.E[i] * eV / (m_e * c * c);
    float v = c * sqrt(1.0 - 1.0 / (gamma * gamma));

    // Calculate the probability of a collision
    float arg = ndens * total_sigma * v * 1e2 * a.dt[i]; // ndens in cm-3, sigma in cm2, v in cm/s
    float p_collision = 1 - exp(-arg); // ndens in cm-3, sigma in cm2, v in cm/s
    // Generate a random number
   float r = getRandom(i);

    // Check if the electron collides
    return r < p_collision;
}

__host__ __device__
int get_collision_type(const SimParams&    p,
                       float*              sigmas,
                       float               total_sigma,
                       int                 tid)
{
    float r = getRandom(tid);

    float cum_prob = 0.0f;
    for (int i = 0; i < p.nCollisions; ++i) {
        cum_prob += sigmas[i] / total_sigma;
        if (r < cum_prob) return i;
    }
    return p.elastic; // Default to elastic collision if no type is found
}

__device__ __host__
void rotate_velocity(float& vz, float& vy, float theta)
{

    float vz_new = vz * cos(theta) - vy * sin(theta);
    float vy_new = vz * sin(theta) + vy * cos(theta);
    vz = vz_new;
    vy = vy_new;
}

__device__ __host__
void isotropic_scattering(float& vz, float& vy, const int i)
{
    // Get the angle
    float r = getRandom(i);
    float theta = r * 2 * pi;

    // Rotate the velocity
    rotate_velocity(vz, vy, theta);
}

__device__ __host__
void screened_rutherford_scattering(const int e,
                                    const SimParams& p, 
                                    DeviceArrays a,
                                    const int i) {

    if (e < 0 || e >= p.nE) return;

    // Sample theta using the Walker alias method
    // Generate a random integer index j in the range [0, nTheta)
    float r = getRandom(i);
    int j = (int)(r * p.nTheta);
    if (j >= p.nTheta) j = p.nTheta - 1;
    // Generate a random number in [0,1]
    float r2 = getRandom(i);
    int idx = j * p.nE + e;
    int t = (r2 < a.prob[idx]) ? j : a.alias[idx];
    // Now sample uniformly within the bin to reduce discretization effects
    // t = p.nTheta - 1 ==> 0 rad (i.e. forward scattering)  so we need to reverse the index
    int k = (p.nTheta - 1 - t);        // k=0 is most forward bin
    float theta_lo = (k    ) * (pi / p.nTheta);
    float theta_hi = (k + 1) * (pi / p.nTheta);
    float theta = theta_lo + getRandom(i) * (theta_hi - theta_lo);

    // Randomize the sign of the angle
    r = getRandom(i);
    if (r < 0.5) {
        theta = -theta;
    }
    
    // Rotate the velocity
    rotate_velocity(a.vz[i], a.vy[i], theta);
}

__device__ __host__
void scale_velocity_deltaE(float& vz, float& vy, float deltaE, float& E)
{
    if (E < deltaE) {
        vz = 0; vy = 0; deltaE = E; E = 0; // electron lost all its energy
    } else {
        float Enew = (E - deltaE);
        E = Enew; // update the energy in eV
        float gamma  = 1.0f + Enew  * eV / (m_e * c * c);
        float v_new  = c * sqrt(1.0f - 1.0f / (gamma * gamma));
        float scale_factor = v_new / sqrt(vz * vz + vy * vy);
        vz *= scale_factor;
        vy *= scale_factor;
    }
#ifdef __CUDA_ARCH__
    atomicAdd(&Elost_primaries, deltaE); // Atomically add deltaE to the global variable
#else
    // Elost_primaries += deltaE; // For host-side execution
#endif
}

// **************************************************************************************
// Collisions
// **************************************************************************************

// Elastic collision
__device__ __host__
void elastic(const int e, 
             const SimParams& p, 
             DeviceArrays a,
             const int i) {

    screened_rutherford_scattering(e, p, a, i);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.elastic * p.nbinsE + e], 1);
#else
    a.colcount[p.elastic * p.nbinsE + e]++;
#endif
}

// ionisation collision
__device__ __host__
void ionisation(const int e,
                const SimParams& p, 
                DeviceArrays a,
                const int i,
                int* Nactive)
{
#ifdef __CUDA_ARCH__
    atomicAdd(a.colcount + p.ionisation * p.nbinsE + e, 1);
#else
    a.colcount[p.ionisation * p.nbinsE + e]++;
#endif

    float E_ion = a.deltaEs[p.ionisation]; // ionisation threshold energy H2 eV.
    if (a.E[i] < E_ion) {
        return;
    }

    // Partition the energy between the scattered and ejected electrons
    float B  = 8.3; // eV the "known function" for H2
    float r = getRandom(i);
    float E_ej = B * tan(r * atan((a.E[i] - E_ion)/(2*B)));
    float E_scat = a.E[i] - E_ion - E_ej;
    if (E_scat <= 1.0f) E_scat = 1.0f;  // minimum energy for the scattered electron

    //Record the altitude of the ionisation event and the energy of the ejected electron
    int k = static_cast<int>(a.z[i] / p.dz);
    int ej = p.EToe(E_ej);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.nion[k * p.nbinsE + ej], 1);
    atomicAdd(&Elost_primaries, E_ion);
    atomicAdd(&Egained_secondaries, E_ej);
#else
    a.nion[k * p.nbinsE + ej]++;
#endif

    // Elastically scatter the incident electron
    int es = p.EToe(E_scat);
    screened_rutherford_scattering(es, p, a, i);

    // Scale the scattered electron velocity to match its new energy
    float gamma = 1.0 + E_scat * eV / (m_e * c * c);
    float v_scat_new = c * sqrt(1.0 - 1.0 / (gamma * gamma));
    float scale_factor = v_scat_new / sqrt(a.vz[i] * a.vz[i] + a.vy[i] * a.vy[i]);
    a.vz[i] *= scale_factor;
    a.vy[i] *= scale_factor;
    a.E[i] = E_scat; // update the energy in eV
    
}

__device__ __host__
void excitation_triplet_a(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_triplet_a], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_triplet_a * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_triplet_a * p.nbinsE + e]++;
#endif
    
}

__device__ __host__
void excitation_triplet_b(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{  
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_triplet_b], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_triplet_b * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_triplet_b * p.nbinsE + e]++;
#endif
}

__device__ __host__
void excitation_triplet_c(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_triplet_c], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_triplet_c * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_triplet_c * p.nbinsE + e]++;
#endif
}

__device__ __host__
void excitation_triplet_e(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{   
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_triplet_e], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_triplet_e * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_triplet_e * p.nbinsE + e]++;
#endif
}

__device__ __host__
void excitation_singlet_B(const int e,
                          const SimParams& p, 
                          DeviceArrays a,
                          const int i)
{
    screened_rutherford_scattering(e, p, a, i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_singlet_B], a.E[i]);
    // Record the altitude and energy of the excitation event
    int k = static_cast<int>(a.z[i] / p.dz);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.nexB[k * p.nbinsE + e], 1);
#else
    a.nexB[k * p.nbinsE + e]++;
#endif
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_singlet_B * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_singlet_B * p.nbinsE + e]++;
#endif
}

__device__ __host__
void excitation_singlet_C(const int e,
                          const SimParams& p, 
                          DeviceArrays a,
                          const int i)
{
    screened_rutherford_scattering(e, p, a, i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_singlet_C], a.E[i]);
    // Record the altitude and energy of the excitation event
    int k = static_cast<int>(a.z[i] / p.dz);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.nexC[k * p.nbinsE + e], 1);
#else
    a.nexC[k * p.nbinsE + e]++;
#endif
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_singlet_C * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_singlet_C * p.nbinsE + e]++;
#endif
}

__device__ __host__
void excitation_singlet_EF(const int e,
                          const SimParams& p, 
                          DeviceArrays a,
                          const int i)
{
    screened_rutherford_scattering(e, p, a, i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.excitation_singlet_EF], a.E[i]);
    // Record the altitude and energy of the excitation event
    int k = static_cast<int>(a.z[i] / p.dz);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.nexEF[k * p.nbinsE + e], 1);
#else
    a.nexEF[k * p.nbinsE + e]++;
#endif
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.excitation_singlet_EF * p.nbinsE + e], 1);
#else
    a.colcount[p.excitation_singlet_C * p.nbinsE + e]++;
#endif
}

__device__ __host__
void vibrational(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{ 
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.vibrational], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.vibrational * p.nbinsE + e], 1);
#else
    a.colcount[p.vibrational * p.nbinsE + e]++; 
#endif
}

__device__ __host__
void rotational(const int e,
                          const SimParams& p, 
                          DeviceArrays a, 
                          const int i)
{
    isotropic_scattering(a.vz[i], a.vy[i], i);
    scale_velocity_deltaE(a.vz[i], a.vy[i], a.deltaEs[p.rotational], a.E[i]);
#ifdef __CUDA_ARCH__
    atomicAdd(&a.colcount[p.rotational * p.nbinsE + e], 1);
#else
    a.colcount[p.rotational * p.nbinsE + e]++;
#endif
}
// **************************************************************************************
// Main kernel
// **************************************************************************************

__device__ __host__
void processElectron(int i,
                        DeviceArrays a, 
                        const SimParams p,
                        ElectronLossCounters* N) {

    // Advance the electron
    a.z[i] += a.vz[i] * a.dt[i];
    a.y[i] += a.vy[i] * a.dt[i];
 
    // Get the electron's energy

    int e = p.EToe(a.E[i]);

    // Check if the electron is still alive. If not, increment the total number of electrons 
    // simulated
    if (a.z[i] < 0 || a.z[i] > p.Z1 || a.E[i] < p.Emin || a.y[i] < -10000e3 || a.y[i] > 10000e3) {
#ifdef __CUDA_ARCH__
        if (a.z[i] < 0.0f) {
            atomicAdd(&N->Nlost_bottom, 1);
        } else if (a.z[i] > p.Z1) {
            atomicAdd(&N->Nlost_top, 1);
        } else if (a.E[i] < p.Emin) {
            atomicAdd(&N->Nlost_energy, 1);
        } else {
            atomicAdd(&N->Nlost_sides, 1);
        }
        // decrement active‐count atomically:
        atomicSub(&N->Nactive, 1);
        // Increment the total energy lost by dead primaries:
        atomicAdd(&Edied_primaries, a.E[i]);
#else
        if (a.z[i] < 0.0f) {
            N->Nlost_bottom++;
        } else if (a.z[i] > p.Z1) {
            N->Nlost_top++;
        } else if (a.E[i] < p.Emin) {
            N->Nlost_energy++;
        } else {
            N->Nlost_sides++;
        }
        N->Nactive--;
#endif
        a.alive[i] = 0;
        return;
    }

    // Get the neutral density at the electron's altitude
    int zinx = static_cast<int>(a.z[i] / p.dz);
    zinx = max(0, min(zinx, p.nbinsZ - 1));
    float ndens = a.nH2[zinx]; // cm^-3

    // Get the cross sections for all collision types and the total for this energy
    float totsig = a.total_sigma[e];

    // Each thread gets its own small array of size p.nCollisions:
    // (small enough to live on the stack)
    float sigmas_local[p.nCollisions];  // e.g. 11 entries

    // Fill the local array from a.sigmasE:
    for (int j = 0; j < p.nCollisions; ++j) {
        sigmas_local[j] = a.sigmasE[e * p.nCollisions + j];
    }

    // Check if the electron collides
    bool collision = test_collision(ndens, totsig, a, i);

    if (collision) {
        // If the electron collides, determine the collision type
        int coltype = get_collision_type(p, sigmas_local, totsig, i);

        switch (coltype) {
            case p.elastic:
                elastic(e, p, a, i);
                break;
            case p.ionisation:
                ionisation(e, p, a, i, &(N->Nactive));
                break;
            case p.excitation_triplet_a:
                excitation_triplet_a(e, p, a, i);
                break;
            case p.excitation_triplet_b:
                excitation_triplet_b(e, p, a, i);
                break;
            case p.excitation_triplet_c:
                excitation_triplet_c(e, p, a, i);
                break;
            case p.excitation_triplet_e:
                excitation_triplet_e(e, p, a, i);
                break;
            case p.excitation_singlet_B:
                excitation_singlet_B(e, p, a, i);
                break;
            case p.excitation_singlet_C:
                excitation_singlet_C(e, p, a, i);
                break;
            case p.vibrational:
                vibrational(e, p, a, i);
                break;
            case p.rotational:
                rotational(e, p, a, i);
                break;
            case p.excitation_singlet_EF:
                excitation_singlet_EF(e, p, a, i);
                break;
            default:
                break;
        }
    }

    // Adapt the timestep for the next iteration based on the electron's speed and the collision probability, to improve efficiency.
    float speed = sqrt(a.vz[i] * a.vz[i] + a.vy[i] * a.vy[i]) * 1e2; // cm/s
    float maxdt = 0.01 / (speed * ndens * totsig); // all cm -ln(0.99)/nsv
    a.dt[i] = (maxdt < p.maxdt) ? maxdt : p.maxdt;

}

__global__
void processElectronsGPU(const DeviceArrays    a,
                         const SimParams       p,
                         ElectronLossCounters* N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p.N) return;
    if (a.alive[i] == 0) {
        return;
    }
     // If here, the electron is still alive
   processElectron(i, a, p, N);
}

__host__ 
void processElectronsCPU(const DeviceArrays    a,
                         const SimParams       p,
                         ElectronLossCounters* N)
{
    for (int i = 0; i < p.N; ++i) {
        if (a.alive[i] == 1) {
            processElectron(i, a, p, N);
        }
    }
}

// Kernel to initialize curandStates
__global__ void setup_kernel(curandState *states, unsigned long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    // Each thread gets the same seed, a different sequence number, no offset
    curand_init(seed,idx, 0, &states[idx]);
}

// **************************************************************************************
// Main function
// **************************************************************************************

void runPrimariesKernel(DeviceArrays d_arrs, DeviceArrays h_arrs, SimParams p) {



    // Allocate memory for the counters on the device
    // This needed to be done on the device because the kernel will modify it
    ElectronLossCounters* d_counters;
    cudaMalloc(&d_counters, sizeof(ElectronLossCounters));
    ElectronLossCounters counters = {p.N, 0, 0, 0, 0};
    cudaMemcpy(d_counters, &counters, sizeof(ElectronLossCounters), cudaMemcpyHostToDevice);

    float Einit = 0.0f;
    for (int i = 0; i < p.N; ++i) {
        Einit += h_arrs.E[i];
    }

    // SimParams p does not need to be copied to the device as it is a small, read-only plain old data structure
    // so can be passed by value to the kernel.


    // Initialize the seed for random number generation
    const int numrands = 1 << 22;  // a way of writing 2^22 == 1048576^2
    curandState *d_states;
    cudaMalloc(&d_states, numrands * sizeof(curandState));

    // Set up the grid and block sizes for one thread per random number
    int threadsPerBlock = 256;
    int blocksPerGrid = (numrands + threadsPerBlock - 1) / threadsPerBlock;
    // --- Set up the cuRAND states on the device ---
    unsigned long seed = static_cast<unsigned long>(time(NULL));
    setup_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_states, seed);
    cudaDeviceSynchronize();

    // Copy the pointer to the global device variable gStates
    cudaMemcpyToSymbol(gStates, &d_states, sizeof(curandState*));

    // Reset device-side energy counters so repeated runs start from zero
    float zero = 0.0f;
    cudaMemcpyToSymbol(Elost_primaries, &zero, sizeof(float));
    cudaMemcpyToSymbol(Edied_primaries, &zero, sizeof(float));
    cudaMemcpyToSymbol(Egained_secondaries, &zero, sizeof(float));
    cudaMemcpyToSymbol(Ealive, &zero, sizeof(float));
    // End of cuRAND setup

    blocksPerGrid = (p.N + threadsPerBlock - 1) / threadsPerBlock;

    // Main loop
    bool advance = true;
    long j = -1;
    std::cout << "Number of active electrons: " << counters.Nactive <<  std::endl;    

    while (advance == true) {
        j++;

        if (j % 10000 == 0) {
            std::cout << "\rNumber of active electrons: " << counters.Nactive << ", " << counters.Nlost_top/(double) (p.N)*100 << "% out the top." << std::flush;
        }

        // Launch the kernel with one thread
        processElectronsGPU<<<blocksPerGrid, threadsPerBlock>>>(d_arrs, p, d_counters);

        // Check launch:
        // CUDA_CHECK(cudaGetLastError());

        // // Check execution:
        // CUDA_CHECK(cudaDeviceSynchronize());

        // Copy the counters back to the host
        cudaMemcpy(&counters, d_counters, sizeof(ElectronLossCounters), cudaMemcpyDeviceToHost);
        
        if (counters.Nactive <= 0) { // If there are no more active electrons, stop the simulation
            std::cout << "\nReached minimum number of active electrons. Stopping...\n";
            advance = false; // Stop the GPU execution
        }

        // if (j > 80000) advance = false; // For testing purposes, stop after 100 iterations
    }
    

    // Copy the final arrays back to the host
    size_t size_E         = p.N * sizeof(float);
    cudaMemcpy(h_arrs.E, d_arrs.E, size_E, cudaMemcpyDeviceToHost);
    size_t size_alive      = p.N * sizeof(int);
    cudaMemcpy(h_arrs.alive, d_arrs.alive, size_alive, cudaMemcpyDeviceToHost);
    size_t size_nion      = p.nbinsZ * p.nbinsE * sizeof(int);
    cudaMemcpy(h_arrs.nion, d_arrs.nion, size_nion, cudaMemcpyDeviceToHost);

    size_t size_theta = p.nTheta * p.nE * sizeof(int);
    cudaMemcpy(h_arrs.theta_sampled, d_arrs.theta_sampled, size_theta, cudaMemcpyDeviceToHost);
    size_t size_colcount = p.nCollisions * p.nbinsE * sizeof(int);
    cudaMemcpy(h_arrs.colcount, d_arrs.colcount, size_colcount, cudaMemcpyDeviceToHost);
    size_t size_nexB = p.nbinsZ * p.nbinsE * sizeof(int);
    cudaMemcpy(h_arrs.nexB, d_arrs.nexB, size_nexB, cudaMemcpyDeviceToHost);
    size_t size_nexC = p.nbinsZ * p.nbinsE * sizeof(int);
    cudaMemcpy(h_arrs.nexC, d_arrs.nexC, size_nexC, cudaMemcpyDeviceToHost);
    size_t size_nexEF = p.nbinsZ * p.nbinsE * sizeof(int);
    cudaMemcpy(h_arrs.nexEF, d_arrs.nexEF, size_nexEF, cudaMemcpyDeviceToHost);
    // Copy the total energy lost
    float Elost_primaries_host, Edied_primaries_host, Egained_secondaries_host;
    cudaMemcpyFromSymbol(&Elost_primaries_host, Elost_primaries, sizeof(float), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&Edied_primaries_host, Edied_primaries, sizeof(float), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&Egained_secondaries_host, Egained_secondaries, sizeof(float), 0, cudaMemcpyDeviceToHost);


    // Free all device memory
    cudaFree(d_counters);
    cudaFree(d_states);
    cudaFree(d_arrs.nH2);
    cudaFree(d_arrs.sigmasE);
    cudaFree(d_arrs.total_sigma);
    cudaFree(d_arrs.sigmas);
    cudaFree(d_arrs.prob);
    cudaFree(d_arrs.alias);
    cudaFree(d_arrs.dt);
    cudaFree(d_arrs.z);
    cudaFree(d_arrs.y);
    cudaFree(d_arrs.vz);
    cudaFree(d_arrs.vy);
    cudaFree(d_arrs.alive);
    cudaFree(d_arrs.dt);
    cudaFree(d_arrs.E);
    // cudaFree(d_arrs.zfinal);
    // cudaFree(d_arrs.yfinal);
    cudaFree(d_arrs.nion);
    // cudaFree(d_arrs.nionz);
    cudaFree(d_arrs.theta_sampled);
    cudaFree(d_arrs.colcount);
    cudaFree(d_arrs.nexB);
    cudaFree(d_arrs.nexC);
    cudaFree(d_arrs.nexEF);

    // Accumulate the energy of alive electrons remaining at the end of the simulation
    double total_alive_energy = 0.0;
    for (int i = 0; i < p.N; ++i) {
        if (h_arrs.alive[i] == 1) {
            total_alive_energy += h_arrs.E[i];
        }
    }

    // Print a report
   cout << "----------------------------------------" << endl;
    cout << "Initial energy in primaries: " << Einit/1e6 << " MeV" << endl;
    cout << "Number of electrons lost at the top boundary: " << counters.Nlost_top << " (" << counters.Nlost_top/(double) (p.N)*100 << "%)" << endl;
    cout << "Number of electrons lost at the bottom boundary: " << counters.Nlost_bottom << " (" << counters.Nlost_bottom/(double) (p.N)*100 << "%)" << endl;
    cout << "Number of electrons lost at the sides: " << counters.Nlost_sides << " (" << counters.Nlost_sides/(double) (p.N)*100 << "%)" << endl;
    cout << "Number of electrons lost due to energy: " << counters.Nlost_energy << " (" << counters.Nlost_energy/(double) (p.N)*100 << "%)" << endl;
    cout << "Total energy lost by primary electrons via collisions: " << Elost_primaries_host/1e6 << " MeV" << endl;
    cout << "Total energy moved from primary to secondary electrons: " << Egained_secondaries_host/1e6 << " MeV" << endl;
    cout << "Total energy removed from primaries: " << (Elost_primaries_host + Egained_secondaries_host)/1e6 << " MeV" << endl;
    cout << "Total energy of primary electrons as they die (e.g. low E): " << Edied_primaries_host/1e6 << " MeV" << endl;
    cout << "Total energy of electrons still alive at end: " << total_alive_energy / 1e6 << " MeV" << endl;
    cout << "Total energy accounted for: " << (Elost_primaries_host + Edied_primaries_host + Egained_secondaries_host + total_alive_energy)/1e6 << " MeV" << endl;
    cout << "----------------------------------------" << endl;

}
