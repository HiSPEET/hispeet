load("CMake")
load("VTK/8.1.0-intel-2018a-Python-3.6.4")
load("impi")
load("VTune")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")

