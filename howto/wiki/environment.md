# Environment Wiki

This wiki explains how to configure the environment for building and running HiSPEET.The simplest approach is to set the environment variables `CC`, `CXX` and `FC`  to the MPI compiler scripts for C, C++ and Fortran, respectively, e.g.

      export CC=mpicc
      export CXX=mpicxx
      export FC=mpifort

However, systems offer a variety of different compilers, MPI implementations and libraries. In this case you need to select specific computing environment witch all ingredients matching to each other. Moreover, you may want to define different environments, one with GCC and  OpenMPI and, maybe, another with Intel compilers and Intel MPI.

TBC