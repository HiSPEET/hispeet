-- VTK with foss/2020b not yet available :(
load("CMake")
load("VTK/8.2.0-foss-2020a-Python-3.8.2")
load("foss/2020b")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

