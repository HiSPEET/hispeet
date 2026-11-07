# Environment Wiki

This wiki explains how to configure the environment for building and running HiSPEET.

## Minimum

1. Set the environment variables `CC`, `CXX` and `FC`  to the MPI compiler scripts for C, C++ and Fortran, e.g.
   
    ```shell
    export CC=mpicc
    export CXX=mpicxx
    export FC=mpifort
    ```
    
2. Take care that the required libraries are installed and can be found
   - BLAS and LAPACK
   - HDF5
   - FFTW3

## HPC sytems at ZIH @ TU Dresden

1. Install the `lua` environment module provided in the _HiSPEET_ repository, e.g.
   
    ```shell
    # in HiSPEET root directory:
    mkdir -p ~/modulefiles/hispeet
    cp -r wiki/environment/barnard ~/modulefiles/hispeet/
    cp -r wiki/environment/capella ~/modulefiles/hispeet/
    ```

2. To make the modules visible, add the following line to your `~/.bashrc`
   
    ```shell
    module use $HOME/modulefiles
    ```

3. To use, e.g., the GNU compilers on barnard, load the corresponding module
   
    ```shell
    module load hispeet/barnard/gnu/14
    ```

## Ubuntu

The same approach is available with Ubuntu, see [README.md](wiki/environment/ubuntu/README.md) in `wiki/environment/ubuntu` 
