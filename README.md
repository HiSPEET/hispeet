# HiSPEET
High-order spectral element techniques.

## About

*HiSPEET* is a framework for developing high-order spectral element methods, with a special focus on Comuputational Fluid Dynamics. HiSPEET builds on unstructured meshes of hexahedral elements with polynomial degrees ranging from 1 to 32 and beyond. It provides fast multigrid solvers and supports adaptive mesh refinement based on the full approximation storage scheme of Achi Brandt. The framework is parallelized using MPI and has been run on thousands of cores.

*HiSPEET* has been developed at the Chair of Fluid Mechanics  at the Insititute of Fluid Mechanics at TU Dresden, and is released under the under the **GNU General Public License v3.0**. For the full license terms see the included [license file](LICENSE.md).

## Dependencies

*HiSPEET* is tested on various Linux distributions and macOS. For installation, the following prerequisites are required: 

- [Git](https://git-scm.com/)
- [CMake](https://cmake.org)
- Fortran and C/C++ Compilers ([GNU](https://gcc.gnu.org/) recommended)
- [BLAS](https://www.netlib.org/blas) and [LAPACK](https://www.netlib.org/lapack)
- [MPI](https://www.mpi-forum.org)
- [HDF5](https://www.hdfgroup.org)
- [FFTW](https://www.fftw.org)

Additionally, [LIBXSMM](https://github.com/libxsmm) is used for the fast evaluation of tensor-product operators on CPUs. This library is installed by *HiSPEET* itself. Results can be visualized using [Python](https://www.python.org) and [Matplotlib](https://matplotlib.org) as well as [ParaView](https://www.paraview.org/). [ford](https://github.com/Fortran-FOSS-Programmers/ford) is utilized for code documentation.

See here how to configure [Ubuntu](./wiki/ubuntu.md) for use with _HiSPEET_. 

## Download

_HiSPEET_ can be downloaded from the GitLab server of TU Chemnitz via

     git clone https://github.com/HiSPEET/hispeet.git

This creates a clone of the git repository in the directory `hispeet`. Next initialize the external libraries which are incorporated as submodules:

      cd hispeet
      git submodule update --init

## Building

For building _HiSPEET_ create a build directory and issue the `cmake` and `make` command from there, e.g:

      mkdir build
      cd build
      cmake ..
      make

Note that all files are generated in the build directory, while the source directories remain unchanged. This allows build the code with different options or even different compilers. For example, you may create a multithreaded version using

      mkdir build-omp
      cd build-omp
      cmake -D OpenMP=1 ..
      make

For further options see the [build](./wiki/build) wiki.
