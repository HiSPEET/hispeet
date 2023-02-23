load("CMake/3.23.1-GCCcore-11.3.0")
load("VTK/9.1.0-foss-2021b")
load("FFTW/3.3.10-gompi-2021b")
load("ParMETIS/4.0.3-gompi-2021b")
load("foss/2021b")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

