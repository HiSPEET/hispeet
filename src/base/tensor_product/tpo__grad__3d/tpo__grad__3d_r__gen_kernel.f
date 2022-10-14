!> summary:  Generic gradient of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Grad_R_Gen_RWP(np, ne, Ms, Ds, dx, u, up, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ms(np)           !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: dx(3)            !< element extensions
  real(RWP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RWP), intent(in)  :: up(np,np,6,ne)   !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

  optional :: Ms, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: g(3), tmp1, tmp2, tmp3
  integer   :: e, f, i, j, k, p, n
  logical   :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. present(up)

  ! metric coefficients
  g = 2 / dx

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! v1 = du/dx1, v2 = du/dx2 .................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + Ds(i,i) * u(p,j,k,e)
        tmp2 = tmp2 + Ds(j,p) * u(i,p,k,e)
      end do
      v(i,j,k,e,1) = g(1) * tmp1
      v(i,j,k,e,2) = g(2) * tmp2

    end do
    end do
    end do

    ! v3 = du/dx3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp3 = 0
      do p = 1, np
        tmp3 = tmp3 + Ds(k,p) * u(i,j,p,e)
      end do
      v(i,j,k,e,3) = g(3) * tmp3

    end do
    end do
    end do

    if (fluxes) then

      ! faces 1 and 2
      i =  1
      n = -1
      do f = 1, 2
        tmp1 = n / (dx(1) * Ms(i))
        do k = 1, np
        do j = 1, np
          v(i,j,k,e,1) = v(i,j,k,e,1) + tmp1 * (up(j,k,f,e) - u(i,j,k,e))
        end do
        end do
        i = np
        n = 1
      end do

      ! faces 3 and 4
      j =  1
      n = -1
      do f = 3, 4
        tmp2 = n / (dx(2) * Ms(j))
        do k = 1, np
        do i = 1, np
          v(i,j,k,e,2) = v(i,j,k,e,2) + tmp2 * (up(i,k,f,e) - u(i,j,k,e))
        end do
        end do
        j = np
        n = 1
      end do

      ! faces 5 and 6
      k =  1
      n = -1
      do f = 5, 6
        tmp3 = n / (dx(3) * Ms(k))
        do j = 1, np
        do i = 1, np
          v(i,j,k,e,3) = v(i,j,k,e,3) + tmp3 * (up(i,j,f,e) - u(i,j,k,e))
        end do
        end do
        k = np
        n = 1
      end do

    end if

  end do

!===============================================================================

end subroutine TPO_Grad_R_Gen_RWP
