load("release/2026")
load("GCC/14.3.0")
load("OpenMPI/5.0.8")
load("HDF5/1.14.6")
load("ParMETIS/4.0.3")
load("CMake/4.0.3")
load("GSL/2.8")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

