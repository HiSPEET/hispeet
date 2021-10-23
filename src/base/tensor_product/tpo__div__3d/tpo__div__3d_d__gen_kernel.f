!> summary:  Generic divergence of a vector field: 3D curvilinear
!> author:   Jerome Michel,Joerg Stiller, Erik Pfister
!> date:     2021/08/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Div_D_Gen_RWP(np, ne, Ds, Ji, u, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: Ji(:,:,:,:,:,:)  !< inverse Jacobi matrix
  real(RWP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
  real(RWP), intent(out) :: v(np,np,np,ne)   !< element-wise divergence of u

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: A(np,np)
  real(RWP) :: r, s, t

  integer :: e, i, j, k, p


  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  !-----------------------------------------------------------------------------
  ! evaluation

  ! v = div(u) = (JiT * grad_xi) * u

  !$omp do private(e)
  do e = 1, ne

      ! direction 1 ............................................................
      do k = 1, np
      do j = 1, np
      do i = 1, np

        r = 0
        s = 0
        t = 0

        ! [r,s,t]^T = [dv1/dxi1 , dv1/dxi2 , dv1/dxi3]
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,1)
          s = s + A(p,j) * u(i,p,k,e,1)
          t = t + A(p,k) * u(i,j,p,e,1)
        end do

        ! v = [JiT_11 , JiT_12 , JiT_13] * [r,s,t]^T
        v(i,j,k,e) =  r * Ji(i,j,k,e,1,1) &
                    + s * Ji(i,j,k,e,2,1) &
                    + t * Ji(i,j,k,e,3,1)

      end do
      end do
      end do


      ! direction 2 ............................................................
      do k = 1, np
      do j = 1, np
      do i = 1, np

        r = 0
        s = 0
        t = 0

        ! [r,s,t]^T = [dv2/dxi1 , dv2/dxi2 , dv2/dxi3]
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,2)
          s = s + A(p,j) * u(i,p,k,e,2)
          t = t + A(p,k) * u(i,j,p,e,2)
        end do

        ! v = v + [JiT_21 , JiT_22 , JiT_23] * [r,s,t]^T
        v(i,j,k,e) =  v(i,j,k,e) + r * Ji(i,j,k,e,1,2) &
                                 + s * Ji(i,j,k,e,2,2) &
                                 + t * Ji(i,j,k,e,3,2)

      end do
      end do
      end do

      ! direction 3 ............................................................
      do k = 1, np
      do j = 1, np
      do i = 1, np

        r = 0
        s = 0
        t = 0

        ! [r,s,t]^T = [dv3/dxi1 , dv3/dxi2 , dv3/dxi3]
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,3)
          s = s + A(p,j) * u(i,p,k,e,3)
          t = t + A(p,k) * u(i,j,p,e,3)
        end do

        ! v = v + [JiT_31 , JiT_32 , JiT_33] * [r,s,t]^T
        v(i,j,k,e) =  v(i,j,k,e) + r * Ji(i,j,k,e,1,3) &
                                 + s * Ji(i,j,k,e,2,3) &
                                 + t * Ji(i,j,k,e,3,3)

      end do
      end do
      end do

  end do
  !$omp end do
!===============================================================================

end subroutine TPO_Div_D_Gen_RWP
