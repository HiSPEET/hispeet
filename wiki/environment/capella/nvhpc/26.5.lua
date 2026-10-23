-- NVIDIA HPC SDK with its MPI (HPC-X OpenMPI) and BLAS/LAPACK
-- * FFTW of the GCC toolchain (C library + fftw3.f03)
-- * no HDF5 for NVHPC in release/2026: own serial build with Fortran + HL
--   (GCC-built HDF5 Fortran modules cannot be read by nvfortran)
-- * CMake 4 needs GCCcore/14.3.0, which NVHPC replaces by 15.2.0
--   -> use the CMake binary directly; the Lmod warning about GCCcore/14.3.0
--   not being loaded is expected and harmless
-- * no ParMETIS
-- * GPU offloading: cmake -D OpenMP=1 -D OFFLOAD=1 [-D GPU_ARCH=cc90]
load("release/2026")
load("GCC/14.3.0")
load("FFTW/3.3.10")
load("NVHPC/26.5-CUDA-13.2.0")

local hdf5 = pathJoin(os.getenv("HOME"), "Numerics/opt/hdf5-1.14.6-nvhpc-26.5")
setenv("HDF5_ROOT", hdf5)
prepend_path("LD_LIBRARY_PATH", pathJoin(hdf5, "lib"))
prepend_path("PATH", "/software/genoa/r2026/CMake/4.0.3-GCCcore-14.3.0/bin")

setenv("BLA_VENDOR", "NVHPC")
setenv("CC",  "mpicc")
setenv("CXX", "mpicxx")
setenv("FC",  "mpifort")
