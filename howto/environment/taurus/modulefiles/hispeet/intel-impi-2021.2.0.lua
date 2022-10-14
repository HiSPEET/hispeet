load("CMake/3.16.4-GCCcore-9.3.0")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("FFTW/3.3.8-intel-2020a")
load("impi/2021.2.0-intel-compilers-2021.2.0")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")


