
local MPIHOME = "/opt/openmpi/GNU-8"

prepend_path( "PATH", MPIHOME .. "/bin" )

setenv( "CC" , MPIHOME .. "/bin/mpicc"   )
setenv( "CXX", MPIHOME .. "/bin/mpicxx"  )
setenv( "FC" , MPIHOME .. "/bin/mpifort" )

setenv( "PFUNIT"  , "/opt/pfunit/GNU-8" )
setenv( "PARMETIS", "/opt/metis/GNU-8"  )
setenv( "EXODUS"  , "/opt/exodus/GNU-8" )
