program Prep

#ifdef __LIBXSMM__
  integer :: LIBXSMM = 1
  print *, '__LIBXSMM__ defined'
  print *, 'LIBXSMM =', LIBXSMM
#else
  print *, '__LIBXSMM__ undefined'
#endif

end program Prep
