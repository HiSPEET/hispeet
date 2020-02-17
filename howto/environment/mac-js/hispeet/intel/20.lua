
local MPIHOME = "/opt/openmpi/INTEL-20"

prepend_path( "PATH", MPIHOME .. "/bin" )

setenv( "CC" , MPIHOME .. "/bin/mpicc"   )
setenv( "CXX", MPIHOME .. "/bin/mpicxx"  )
setenv( "FC" , MPIHOME .. "/bin/mpifort" )

setenv( "PFUNIT"  , "/opt/pfunit/INTEL-20" )
setenv( "PARMETIS", "/opt/metis/INTEL-20"  )
setenv( "EXODUS"  , "/opt/exodus/INTEL-20" )
