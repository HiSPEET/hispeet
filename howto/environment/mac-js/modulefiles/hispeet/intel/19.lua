
local MPIHOME = "/opt/openmpi/INTEL-19"

prepend_path( "PATH", MPIHOME .. "/bin" )

setenv( "CC" , MPIHOME .. "/bin/mpicc"   )
setenv( "CXX", MPIHOME .. "/bin/mpicxx"  )
setenv( "FC" , MPIHOME .. "/bin/mpifort" )

setenv( "PFUNIT"  , "/opt/pfunit/INTEL-19" )
setenv( "PARMETIS", "/opt/metis/INTEL-19"  )
setenv( "EXODUS"  , "/opt/exodus/INTEL-19" )
