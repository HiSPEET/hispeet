load("CMake")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("impi")
load("VTune")

setenv("ENV_MODULE_COMPILER", "hispeet/intel-impi")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")

