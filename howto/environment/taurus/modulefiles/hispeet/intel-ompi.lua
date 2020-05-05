load("CMake/3.13.3-GCCcore-8.2.0")
load("VTK/8.1.0-intel-2018a-Python-3.6.4")
load("OpenMPI/3.1.3-iccifort-2019.1.144-GCC-8.2.0-2.31.1")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

