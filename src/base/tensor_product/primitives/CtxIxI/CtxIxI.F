#ifdef __LIBXSMM__

#  ifdef __INTEL_COMPILER
#    include "CtxIxI_nc02-07_intel.F"
#  else
#    include "CtxIxI_nc02-07_gnu.F"
#  endif

#else

#  ifdef __INTEL_COMPILER
#    include "CtxIxI_nc02-33_intel.F"
#  else
#    include "CtxIxI_nc02-33_gnu.F"
#  endif

#endif

