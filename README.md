# HiSPEET Primer
This document explains how to make your first steps with _HiSPEET_.

## Prerequisites

In order to install _HiSPEET_ you need a Unix-like operating system with

- bash, git, cmake, make
- Fortran, C and C++ compilers, Python
- MPI and VTK libraries

and, ideally,

- Paraview for visualization
- FFTW for some postprocessing tools.

See here how to set up [Git](./howto/wiki/git.md) and how to configure [Ubuntu](./howto/wiki/ubuntu.md) or macOS for use with _HiSPEET_.

## Download

*needs to be revised*
_HiSPEET_ can be obtained via [anonymous access](https://fusionforge.zih.tu-dresden.de/ scm/?group_id=747)  from the ZIH FusionForge server

     git clone https://scm.fusionforge.zih.tu-dresden.de/anonscm/git/hispeet/hispeet.git

This creates a clone of the git repository in the directory `hispeet`.  Next initialize the external libraries which are incorporated as submodules

      cd hispeet
      git submodule update --init

## Environment

Building _HiSPEET_ requires a suitable computing environment. Often it is sufficient to set the environment variables `CC`, `CXX` and `FC` to the MPI compiler scripts for C, C++ and Fortran, respectively, e.g.

      export CC=mpicc
      export CXX=mpicxx
      export FC=mpifort

If your system offers different sets of compilers and MPI implementations, it is advisable to configure the environment for choosing the right one. See [environment](./howto/wiki/environment) wiki for details and examples.

## Building

For building _HiSPEET_ create a build directory and issue the `cmake` and `make` command from there

      mkdir build
      cd build
      cmake ..
      make

Note that all files are generated in the build directory, while the source directories remain unchanged. This allows build the code with different options or even different compilers. For example, you may create a multithreaded version using

      mkdir build-omp
      cd build-omp
      cmake -D OpenMP=1 ..
      make

For further options see the [build](./howto/wiki/build) wiki.

## Running

TBD

