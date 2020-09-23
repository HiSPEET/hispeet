load("CMake")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("impi")
--load("iimpi/2020a")
load("iccifort/2020.2.254")
load("VTune")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")

