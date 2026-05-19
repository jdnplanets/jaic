# Jupiter Auroral Ionosphere Code (JAIC)

JAIC is an open-source code for modelling electron precipitation and related ionospheric / auroral processes in Jupiter’s upper atmosphere. The hybrid model follows the primary electron population using a Monte Carlo code that runs on an NVIDIA GPU, and computes the contribution of the secondaries using a two-stream approximation.  The model includes modules that compute high resolution far-ultraviolet H2 spectra, the H3+ density using simple ion chemistry, and the resulting Pedersen conductivity and H3+ radiance. 

> **Status:** research code (active development). Interfaces, file formats, and outputs may change.

---

## Release / archive

This repository is archived on **Zenodo**:

- **Zenodo record:** 10.5281/zenodo.18625226

---

## Quick start

### 1) Get the code
If you downloaded a Zenodo ZIP, unpack it and `cd` into the extracted folder.

If you are cloning from a mirror:
```bash
git clone https://github.com/jdnplanets/jaic.git
cd <repo>


### 2) Add the root directory to your environment and add the executable directory to your path, e.g. in a .bashrc file
export JAIC_ROOT=/absolute/path/to/jaic
PATH="$JAIC_ROOT/.:$JAIC_ROOT/bin/.:$PATH"

### 3) Build using the makefile. Note running Precip requires a CUDA/GPU build, otherwise use a CPU build

CPU build:

make CPU=1


CUDA/GPU build:
make

If you swap between these two you will need to include the -B argument to make, to ensure everything is recompiled, e.g.

make -B

Debug build:

make -B DEBUG=1
make -B CPU=1 DEBUG=1

cuda-gdb jaic

### 4) Requirements
CPU build:
C++17 compiler (e.g. g++ or clang)
make

GPU build:
The above plus
NVIDIA GPU
CUDA toolkit (nvcc) and compatible driver
Note: If you are compiling on a cluster, load the appropriate compiler/CUDA modules first.


### 4) Running JAIC
The main executable is:
bin/jaic

If you have added the bin/ directory to your path the code is run by typing

jaic

An example main.cpp file is in the src directory with several sources indicated.
Precip needs to be run once for a given source, yielding results per primary electron
Subsequent calculations require a flux F to be provided, e.g. as shown.
Later modules can either take the results directly from Precip or read in the saved data file if Precip doesn't need to be run again
Only Precip needs the GPU.

## Outputs

Outputs are written under:

${JAIC_ROOT}/out/<module>/<runid>/
where:
<module> refers to the JAIC module <precip>, <h2spec> etc.
<runid> is defined by .runid in the source parameters

The output columns and quantity units are given in the metadata at the top of each file


## Directory layout

include/        Public headers
src/            Source code
bin/            Executables
obj/            Shared objects and dependancy files
data/           Input data tables
out/            Output directory
input/          User-supplied electron spectra
jobs/           HPC submission scripts
logs/           Running logs if necessary with HPC runs


##  User-supplied spectrum
Files should be placed in the input/ directory and consist of two columns (no header).
Column 1: Energy [keV]
Column 2: Differential flux [cm^-2 s^-1 sr^-1 keV^-1]


##  Data
JAIC expects requisite data files under the data/ tree (e.g. precipitation inputs, cross sections, atmospheric profiles). These are included in this repo and details are as follows.

The H2 line emission probabilities are from the MOLAT database:
https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2can94.php
and for further details please see:
Abgrall, H., Roueff, E., Launay, F., Roncin, J.-Y., 1994, Can. J. Phys., 72, 856-865.

The H2 continuum probabilities are also from the MOLAT database:
https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2cont97q5.php
and for further details please see:
Abgrall, H., Roueff, E., Liu, X., Shemansky D.E. 1997, Astroph. J., 481, 557-566

The X->EF Franck-Condon factors were kindly supplied by E. Roueff.

Hydrocarbon photon absorption cross sections are interpolated from Table 5 in:
Parkinson, C. D., J. C. McConnell, L. Ben Jaffel, A. Y. T. Lee, Y. L. Yung, and E. Griffioen. ‘Deuterium Chemistry and Airglow in the Jovian Thermosphere’. Icarus 183, no. 2 (2006): 451–70. https://doi.org/10.1016/j.icarus.2005.09.022.

Hydrocarbon mixing ratios are read from Figure 3b in:
Gladstone, G. Randall, Mark Allen, and Y. L. Yung. ‘Hydrocarbon Photochemistry in the Upper Atmosphere of Jupiter’. Icarus 119, no. 1 (1996): 1–52. https://doi.org/10.1006/icar.1996.0001.

The H2 electron cross sections are read from Figure 12a in:
Hiraki, Y., and C. Tao. ‘Parameterization of Ionization Rate by Auroral Electron Precipitation in Jupiter’. Annales Geophysicae 26, no. 1 (2008): 77–86. https://doi.org/10.5194/angeo-26-77-2008.

The H3+ line list was obtained from the h3ppy website:
https://h3ppy.readthedocs.io/en/master/index.html
originally from the ExoMol database:
https://www.exomol.com/data/molecules/H3_p/1H3_p/NMT/
and 
Neale, Liesl, Steven Miller, and Jonathan Tennyson. ‘Spectroscopic Properties of the H 3 + Molecule: A New Calculated Line List’. The Astrophysical Journal 464 (June 1996): 516. https://doi.org/10.1086/177341.


##  Citation
If you use JAIC in your work, please cite 
Nichols, J.D., Jupiter's auroral ionosphere: Hybrid Monte Carlo, auroral spectrum and conductivity modeling, J. Geophys. Res. (2026)


##  License
MIT License. See LICENSE.txt