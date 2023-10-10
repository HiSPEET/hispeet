
    integer :: i, j, l, m

    m = size(vn,1)
    l = m + 1

    select case(f)
    case(1,2)
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 31_IXS, 51_IXS)
        ve = vn
      case(16_IXS, 35_IXS, 56_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(15_IXS, 34_IXS, 54_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(13_IXS, 32_IXS, 53_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(42_IXS, 61_IXS, 21_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(46_IXS, 65_IXS, 26_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(45_IXS, 64_IXS, 24_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case(43_IXS, 62_IXS, 23_IXS)
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    case(3,4)
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 16_IXS, 24_IXS)
        ve = vn
      case(62_IXS, 56_IXS, 64_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(42_IXS, 46_IXS, 54_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(32_IXS, 26_IXS, 34_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(15_IXS, 13_IXS, 21_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(65_IXS, 53_IXS, 61_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(45_IXS, 43_IXS, 51_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case(35_IXS, 23_IXS, 31_IXS)
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    case default ! 5,6
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 13_IXS, 23_IXS)
        ve = vn
      case(51_IXS, 61_IXS, 62_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(45_IXS, 46_IXS, 56_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(24_IXS, 34_IXS, 35_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(15_IXS, 16_IXS, 26_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(54_IXS, 64_IXS, 65_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(42_IXS, 43_IXS, 53_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case(21_IXS, 31_IXS, 32_IXS)
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    end select
