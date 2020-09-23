load("CMake/3.13.3-GCCcore-8.2.0")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("OpenMPI/3.1.3-iccifort-2019.1.144-GCC-8.2.0-2.31.1")
setenv("ENV_MODULE_COMPILER", "hispeet/intel-ompi")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

