!> summary:  Generic divergence of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Div_R_Gen_RWP(np, ne, Ms, Ds, dx, u, up, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ms(np)           !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: dx(3)            !< element extensions
  real(RWP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
  real(RWP), intent(in)  :: up(np,np,6,ne,3) !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(np,np,np,ne)   !< divergence of u

  optional :: Ms, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: g(3), tmp

  integer :: e, f, i, j, k, p, n
  logical :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. present(up)

  ! metric coefficients
  g = 2 / dx

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! v = du1/dx1 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(i,p) * u(p,j,k,e,1)
      end do
      v(i,j,k,e) = g(1) * tmp

    end do
    end do
    end do

    ! v+= du2/dx2 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(j,p) * u(i,p,k,e,2)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp

    end do
    end do
    end do

    ! v+= du3/dx3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(k,p) * u(i,j,p,e,3)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(3) * tmp

    end do
    end do
    end do

    if (fluxes) then

      ! faces 1 and 2
      i =  1
      n = -1
      do f = 1, 2
        tmp = n / (dx(1) * Ms(i))
        do k = 1, np
        do j = 1, np
          v(i,j,k,e) = v(i,j,k,e) + tmp * (up(j,k,f,e,1) - u(i,j,k,e,1))
        end do
        end do
        i = np
        n = 1
      end do

      ! faces 3 and 4
      j =  1
      n = -1
      do f = 3, 4
        tmp = n / (dx(2) * Ms(j))
        do k = 1, np
        do i = 1, np
          v(i,j,k,e) = v(i,j,k,e) + tmp * (up(i,k,f,e,2) - u(i,j,k,e,2))
        end do
        end do
        j = np
        n = 1
      end do

      ! faces 5 and 6
      k =  1
      n = -1
      do f = 5, 6
        tmp = n / (dx(3) * Ms(k))
        do j = 1, np
        do i = 1, np
          v(i,j,k,e) = v(i,j,k,e) + tmp * (up(i,j,f,e,3) - u(i,j,k,e,3))
        end do
        end do
        k = np
        n = 1
      end do

    end if

  end do

!===============================================================================

end subroutine TPO_Div_R_Gen_RWP
