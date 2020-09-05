-- REMARK
-- Normally we shold just use `ml foss` to get GCC including OpenMPI, LAPACK etc.

load("CMake")
load("VTK/8.2.0-foss-2020a-Python-3.8.2")
load("OpenMPI/4.0.4")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

-- workaround to ensure that libgfortran.so.4 is found
append_path("LD_LIBRARY_PATH","/sw/installed/GCCcore/7.3.0/lib64")

