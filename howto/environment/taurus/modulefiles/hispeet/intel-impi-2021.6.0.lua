load("CMake/3.23.1-GCCcore-11.3.0")
load("VTK/8.2.0-intel-2020a-Python-3.8.2")
load("FFTW/3.3.10-iimpi-2021b")
load("ParMETIS/4.0.3-iimpi-2020b")
load("impi/2021.6.0-intel-compilers-2022.1.0")

setenv("CC",  "mpicc")
setenv("CXX", "mpiicpc")
setenv("FC",  "mpiifort")

-- setenv( "ParMETIS_ROOT", "?"  )

-- circumvent bug "mpi/pmi2: value not properly terminated in client request"
setenv("I_MPI_PMI_LIBRARY", "/usr/lib64/libpmi2.so")


