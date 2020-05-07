#ifdef __LIBXSMM__

#  ifdef __INTEL_COMPILER
#    include "IxBtxI_na02-07_intel.F"
#  else
#    include "IxBtxI_na02-07_gnu.F"
#  endif

#else

#  ifdef __INTEL_COMPILER
#    include "IxBtxI_na02-33_intel.F"
#  else
#    include "IxBtxI_na02-33_gnu.F"
#  endif

#endif

