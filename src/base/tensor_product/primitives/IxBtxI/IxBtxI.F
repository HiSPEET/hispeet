#ifdef __LIBXSMM__

#  ifdef __INTEL_COMPILER
#    include "IxBtxI_nb02-07_intel.F"
#  else
#    include "IxBtxI_nb02-07_gnu.F"
#  endif

#else

#  ifdef __INTEL_COMPILER
#    include "IxBtxI_nb02-33_intel.F"
#  else
#    include "IxBtxI_nb02-33_gnu.F"
#  endif

#endif

