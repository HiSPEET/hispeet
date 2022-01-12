load("CMake")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("FFTW/3.3.8-intel-2020a")
load("impi/2019.7.217-iccifort-2020.1.217")
load("iccifort/2020.2.254")
load("VTune")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")


