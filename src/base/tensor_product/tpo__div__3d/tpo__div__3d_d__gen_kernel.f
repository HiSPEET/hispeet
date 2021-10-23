!> summary:  Generic divergence of a vector field: 3D curvilinear
!> author:   Jerome Michel, Joerg Stiller
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

  ! v = ∇⋅u = J⁻ᵀ ⋅ (∇ˢ ⋅ u)

  !$omp do private(e)
  do e = 1, ne

      ! direction 1 ............................................................

      do k = 1, np
      do j = 1, np
      do i = 1, np

        r = 0
        s = 0
        t = 0

        ! [r,s,t]ᵀ = ∇ˢu₁
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,1)
          s = s + A(p,j) * u(i,p,k,e,1)
          t = t + A(p,k) * u(i,j,p,e,1)
        end do

        ! v = [J⁻ᵀ₁₁ , J⁻ᵀ₁₂ , J⁻ᵀ₁₃] ⋅ [r,s,t]^T
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

        ! [r,s,t]ᵀ = ∇ˢu₂
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,2)
          s = s + A(p,j) * u(i,p,k,e,2)
          t = t + A(p,k) * u(i,j,p,e,2)
        end do

        ! v += [J⁻ᵀ₂₁ , J⁻ᵀ₂₂ , J⁻ᵀ₂₃] ⋅ [r,s,t]^T
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

        ! [r,s,t]ᵀ = ∇ˢu₃
        do p = 1, np
          r = r + A(p,i) * u(p,j,k,e,3)
          s = s + A(p,j) * u(i,p,k,e,3)
          t = t + A(p,k) * u(i,j,p,e,3)
        end do

        ! v += [J⁻ᵀ₃₁ , J⁻ᵀ₃₂ , J⁻ᵀ₃₃] ⋅ [r,s,t]^T
        v(i,j,k,e) =  v(i,j,k,e) + r * Ji(i,j,k,e,1,3) &
                                 + s * Ji(i,j,k,e,2,3) &
                                 + t * Ji(i,j,k,e,3,3)

      end do
      end do
      end do

  end do
  !$omp end do

end subroutine TPO_Div_D_Gen_RWP
