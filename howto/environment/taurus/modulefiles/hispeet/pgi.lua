load("CMake")
load("VTK/8.1.1-foss-2018b-Python-3.6.6")
load("OpenMPI/4.0.1-PGI-19.4-GCC-8.2.0-2.31.1")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

