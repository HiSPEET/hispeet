local MPIHOME = "/opt/openmpi/GNU-11"

prepend_path( "PATH", MPIHOME .. "/bin" )

-- OpenMPI installed myself
setenv( "CC"    , MPIHOME .. "/bin/mpicc"   )
setenv( "CXX"   , MPIHOME .. "/bin/mpicxx"  )
setenv( "FC"    , MPIHOME .. "/bin/mpifort" )
setenv( "MPIRUN", MPIHOME .. "/bin/mpirun"  )

setenv( "ParMETIS_ROOT", "/opt/metis/GNU-11"  )
--setenv( "Exodus_ROOT"  , "/opt/exodus/GNU-10" )
--setenv( "pFUnit_ROOT"  , "/opt/pfunit/GNU-10" )

