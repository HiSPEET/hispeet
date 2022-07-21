    integer :: i, j, l, m

    m = size(ve,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j,   i)
      end select
    end if
