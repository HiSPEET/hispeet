#ifdef LIBXSMM

#  ifdef __INTEL_COMPILER
#    include "IxQtxI_nq02-07_intel.F"
#  else
#    include "IxQtxI_nq02-07_gnu.F"
#  endif

#else

#  ifdef __INTEL_COMPILER
#    include "IxQtxI_nq02-33_intel.F"
#  else
#    include "IxQtxI_nq02-33_gnu.F"
#  endif

#endif

