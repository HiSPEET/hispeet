    integer :: i, j, l, m

    m = size(vm,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j,   i) = vm(i,j)
      end select
    end if
