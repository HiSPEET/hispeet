-- MPI --

local MPI_HOME = "/opt/homebrew"

setenv( "CC"    , MPI_HOME .. "/bin/mpicc"   )
setenv( "CXX"   , MPI_HOME .. "/bin/mpicxx"  )
setenv( "FC"    , MPI_HOME .. "/bin/mpifort" )
setenv( "MPIRUN", MPI_HOME .. "/bin/mpirun"  )

-- Metis and ParMETIS --

setenv( "ParMETIS_ROOT", "/opt/parmetis/brew-gnu-16" )
