-- REMARK
-- Normally we shold just use `ml foss` to get GCC including OpenMPI, LAPACK etc.

load("CMake")
load("VTK/8.1.1-foss-2018b-Python-3.6.6")
load("Score-P/6.0-gompi-2019b")

setenv("CC",  "scorep-mpicc")
setenv("CXX", "scorep-mpicxx")
setenv("FC",  "scorep-mpif90")

-- workaround to ensure that libgfortran.so.4 is found
append_path("LD_LIBRARY_PATH","/sw/installed/GCCcore/7.3.0/lib64")

