#ifdef LIBXSMM

#  ifdef __INTEL_COMPILER
#    include 'QtxIxI_nq02-07_intel.F'
#  elif
#    include 'QtxIxI_nq02-07_gnu.F'
#  endif

#elif

#  ifdef __INTEL_COMPILER
#    include 'QtxIxI_nq02-33_intel.F'
#  elif
#    include 'QtxIxI_nq02-33_gnu.F'
#  endif

#endif

