load("intel/2023b")
load("iimpi/2023b")
load("HDF5/1.14.3-serial")
load("ParMETIS/4.0.3")
load("VTK/9.3.0")
load("FFTW/3.3.10")

setenv("CC",  "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")

