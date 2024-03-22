load("GCC/11.3.0")
load("OpenMPI/4.1.4")
load("ParMETIS/4.0.3")
load("CMake/3.24.3")
load("VTK/9.2.2")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

