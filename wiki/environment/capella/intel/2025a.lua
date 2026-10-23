-- no ParMETIS for the Intel toolchain in release/2026
load("release/2026")
load("intel/2025a")
load("HDF5/1.14.6")
load("FFTW/3.3.10")
load("CMake/3.31.3")

setenv("CC",  "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC",  "mpiifx")
