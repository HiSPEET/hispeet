
### Prerequisites

     # load SCS5, best put in ~/.bashrc_all
     ml modenv/scs5

     # if conflicting modules are loaded
     module purge


### Intel

     ml CMake
     ml VTK/8.1.0-intel-2018a-Python-3.6.4

     # for use with OpenMPI
     ml OpenMPI/2.1.2-iccifort-2018.1.163-GCC-6.4.0-2.28
     export FC=mpifort
     export CC=mpicc
     export CXX=mpicxx

     # for use with Intel MPI
     ml intel
     export FC=mpiiifort
     export CC=mpiicc
     export CXX=mpiicxx



### GNU

     ml CMake
     ml VTK/8.1.1-foss-2018b-Python-2.7.15
     ml Python/3.6.6-foss-2018b
     ml OpenMPI/3.1.2-GCC-8.2.0-2.30

     export FC=mpifort
     export CC=mpicc
     export CXX=mpicxx

