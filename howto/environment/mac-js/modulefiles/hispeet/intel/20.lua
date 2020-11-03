
local MPIHOME = "/opt/openmpi/INTEL-20"

prepend_path( "PATH", MPIHOME .. "/bin" )

setenv( "CC" , MPIHOME .. "/bin/mpicc"   )
setenv( "CXX", MPIHOME .. "/bin/mpicxx"  )
setenv( "FC" , MPIHOME .. "/bin/mpifort" )

setenv( "FFTW_ROOT"    , "/opt/fftw/INTEL-20"   )
setenv( "Exodus_ROOT"  , "/opt/exodus/INTEL-20" )
setenv( "ParMETIS_ROOT", "/opt/metis/INTEL-20"  )
setenv( "pFUnit_ROOT"  , "/opt/pfunit/INTEL-20" )
