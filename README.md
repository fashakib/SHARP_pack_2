# SHARP pack software package
**SHARP pack**, short for _**Surface Hopping And Ring Polymer package**_, is a highly-parallelized software with a modular structure that is developed to (1) help the development of approximate quantum dynamics methods and (2) simulate non-adiabatic dynamics in condensed phases subject to nuclear quantum effects.

The SHARP Pack documentation is also available on [https://sharppack.readthedocs.io/en/latest/](https://sharppack.readthedocs.io/en/latest/)

## Installation
The SHARP pack software package can be downloaded from the Github repository. Create a folder and git clone this repository.
```
$ git clone https://github.com/fashakib/SHARP_pack.git
```

Once the package has been downloaded, navigate to the directory containing the package and run the following command:
```
$ make
```

## Prerequisite
The SHARP pack is tested under a Linux environment with an Intel Fortran compiler. Path for BLAS, LAPACK, and FFTW3 libraries need to be provided.

## Usage
The SHARP pack software package includes several modules that can be used to perform different types of simulations and analyses. See the SHARP Pack manual for more detail (coming soon).

## Examples
The SHARP pack software package includes several examples that demonstrate how to use the different method(s)/model(s). These examples can be found in the example/ directory of the package.

## Running Simulation
After compiling, **sharp.x** can be called upon to run the simulation in a serial mode
(**ncore 1**). An input file, _**param.in**_, is required to run any SHARP pack simulation.
Alternatively, the SHARP pack can be run using the job submission bash script (see
_**job_script_hpc.sh**_ script in **utility/** directory). Based on **ncore** in _**param.in**_, the code will run as _serial_ on a single core or _parallel_ (trajectories are parallelized in this case) on nodes on the cluster.

## Changes in version 2.0 compared to version 1.0 of the software• Extending the methodology to bead-approximation ring polymer surface hopping (RPSH-BA).• Adding diagonal Born-Oppenheimer Correction (DBOC) implementation.• Option for decoherence correction.• Random seed handling updated for reproducibility.• Added PILE thermostat for PIMD sampling simulations.• New PIMD and runtime outputs: bead- and centroid-level system temperatures.• New option for constraining the centroid to dividing surface (R=0) for initial sampling, shifting the whole ring polymeraccordingly.• Optput histogram distribution of initial sampling positions and momenta.• Extented model system to spin-boson and dboc1 models.• Added Debye/Ohmic spectral density options for spin-boson model.
## Bugs Fixed • Bugs associated with the initial bead position sampling from de Broglie length for spin-boson model are fixed.• Bugs associated with the initial force calculation at t=0 carried out before running dynamics are fixed.
## DocumentationOnline documentation is launched (docs in Sphinx via Readthedocs)

## Contact
For any queries and feedback please contact Dr. Limbu (dil.k.limbu@njit.edu) or Dr. Shakib (shakib@njit.edu).
