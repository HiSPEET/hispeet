-- REMARK
-- Normally we shold just use `ml foss` to get GCC including OpenMPI, LAPACK etc.
-- However, foss/2018b is based one GCC 7, whereas we need GCC 8.

load("CMake")
load("VTK/8.1.1-foss-2018b-Python-3.6.6")
load("OpenMPI/4.0.1-GCC-9.1.0-2.32")

setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")

-- workaround to ensure that libgfortran.so.4 is found
append_path("LD_LIBRARY_PATH","/sw/installed/GCCcore/7.3.0/lib64")

