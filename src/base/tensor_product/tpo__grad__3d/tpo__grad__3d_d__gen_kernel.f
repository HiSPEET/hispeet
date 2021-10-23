!> summary:  Generic gradient of a vector field: 3D Cartesian curvilinear
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Grad_D_Gen_RWP(np, ne, Ds, Ji, u, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: Ji(:,:,:,:,:,:)  !< inverse Jacobi matrix
  real(RWP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RWP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: A(np,np)
  real(RWP) :: r, s, t
  integer   :: e, i, j, k, p

  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    do k = 1, np
    do j = 1, np
    do i = 1, np

      ! [r,s,t]ᵀ = ∇ˢu .........................................................

      r = 0
      s = 0
      t = 0

      do p = 1, np
        r = r + A(p,i) * u(p,j,k,e)
        s = s + A(p,j) * u(i,p,k,e)
        t = t + A(p,k) * u(i,j,p,e)
      end do

      ! v = J⁻ᵀ ⋅ [r,s,t]^T ....................................................

      v(i,j,k,e,1) = r * Ji(i,j,k,e,1,1) &
                   + s * Ji(i,j,k,e,2,1) &
                   + t * Ji(i,j,k,e,3,1)

      v(i,j,k,e,2) = r * Ji(i,j,k,e,1,2) &
                   + s * Ji(i,j,k,e,2,2) &
                   + t * Ji(i,j,k,e,3,2)

      v(i,j,k,e,3) = r * Ji(i,j,k,e,1,3) &
                   + s * Ji(i,j,k,e,2,3) &
                   + t * Ji(i,j,k,e,3,3)

    end do
    end do
    end do

  end do
  !$omp end do

  !-----------------------------------------------------------------------------

end subroutine TPO_Grad_D_Gen_RWP
