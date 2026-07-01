load("release/23.10")
load("GCC/12.2.0")
load("OpenMPI/4.1.4")
load("HDF5/1.14.0")
load("ParMETIS/4.0.3")
load("CMake/3.24.3")
load("VTK/9.2.6")
load("GSL/2.7")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

