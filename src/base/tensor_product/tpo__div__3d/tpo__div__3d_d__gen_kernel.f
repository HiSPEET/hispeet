!> summary:  Generic divergence of a vector field: 3D curvilinear
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Div_D_Gen_RWP(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
  integer,   intent(in)  :: np                  !< num points per direction
  integer,   intent(in)  :: ne                  !< num elements
  real(RWP), intent(in)  :: Ms(np)              !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(np,np)           !< 1D standard diff matrix
  real(RWP), intent(in)  :: Jd(np,np,np,ne)     !< element Jacobian determinant
  real(RWP), intent(in)  :: Ji(np,np,np,ne,3,3) !< inverse Jacobi matrix
  real(RWP), intent(in)  :: a(np,np,6,ne)       !< face area coefficients
  real(RWP), intent(in)  :: n(np,np,6,ne,3)     !< face unit normal vectors
  real(RWP), intent(in)  :: u(np,np,np,ne,3)    !< 3D vector field
  real(RWP), intent(in)  :: up(np,np,6,ne,3)    !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(np,np,np,ne)      !< divergence of u

  optional :: Ms, Jd, a, n, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: c, r, s, t
  integer   :: e, f, i, j, k, p
  logical   :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. &
           present(Jd) .and. &
           present(a)  .and. &
           present(n)  .and. &
           present(up)

  !-----------------------------------------------------------------------------
  ! evaluation

  ! v = ∇⋅u = J⁻ᵀ ⋅ (∇ˢ ⋅ u)

  !$omp do private(e)
  do e = 1, ne

    ! direction 1 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      r = 0
      s = 0
      t = 0

      ! [r,s,t]ᵀ = ∇ˢu₁
      do p = 1, np
        r = r + Ds(i,p) * u(p,j,k,e,1)
        s = s + Ds(j,p) * u(i,p,k,e,1)
        t = t + Ds(k,p) * u(i,j,p,e,1)
      end do

      ! v = [J⁻ᵀ₁₁ , J⁻ᵀ₁₂ , J⁻ᵀ₁₃] ⋅ [r,s,t]^T
      v(i,j,k,e) =  r * Ji(i,j,k,e,1,1) &
                  + s * Ji(i,j,k,e,2,1) &
                  + t * Ji(i,j,k,e,3,1)

    end do
    end do
    end do

    ! direction 2 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      r = 0
      s = 0
      t = 0

      ! [r,s,t]ᵀ = ∇ˢu₂
      do p = 1, np
        r = r + Ds(i,p) * u(p,j,k,e,2)
        s = s + Ds(j,p) * u(i,p,k,e,2)
        t = t + Ds(k,p) * u(i,j,p,e,2)
      end do

      ! v += [J⁻ᵀ₂₁ , J⁻ᵀ₂₂ , J⁻ᵀ₂₃] ⋅ [r,s,t]^T
      v(i,j,k,e) =  v(i,j,k,e) + r * Ji(i,j,k,e,1,2) &
                               + s * Ji(i,j,k,e,2,2) &
                               + t * Ji(i,j,k,e,3,2)

    end do
    end do
    end do

    ! direction 3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np

      r = 0
      s = 0
      t = 0

      ! [r,s,t]ᵀ = ∇ˢu₃
      do p = 1, np
        r = r + Ds(i,p) * u(p,j,k,e,3)
        s = s + Ds(j,p) * u(i,p,k,e,3)
        t = t + Ds(k,p) * u(i,j,p,e,3)
      end do

      ! v += [J⁻ᵀ₃₁ , J⁻ᵀ₃₂ , J⁻ᵀ₃₃] ⋅ [r,s,t]^T
      v(i,j,k,e) =  v(i,j,k,e) + r * Ji(i,j,k,e,1,3) &
                               + s * Ji(i,j,k,e,2,3) &
                               + t * Ji(i,j,k,e,3,3)

    end do
    end do
    end do

    ! fluxes ...................................................................

    if (fluxes) then

      ! faces 1 and 2
      do f = 1, 2
        i = 1 + (np - 1) * (f - 1)
        do k = 1, np
        do j = 1, np
          c = a(j,k,f,e) / (2 * Ms(i ) * Jd(i,j,k,e))
          v(i,j,k,e) = v(i,j,k,e)                                         &
                     + c * n(j,k,f,e,1) * (up(j,k,f,e,1) - u(i,j,k,e,1))  &
                     + c * n(j,k,f,e,2) * (up(j,k,f,e,2) - u(i,j,k,e,2))  &
                     + c * n(j,k,f,e,3) * (up(j,k,f,e,3) - u(i,j,k,e,3))
        end do
        end do
      end do

      ! faces 3 and 4
      do f = 3, 4
        j = 1 + (np - 1) * (f - 3)
        do k = 1, np
        do i = 1, np
          c = a(i,k,f,e) / (2 * Ms(j) * Jd(i,j,k,e))
          v(i,j,k,e) = v(i,j,k,e)                                         &
                     + c * n(i,k,f,e,1) * (up(i,k,f,e,1) - u(i,j,k,e,1))  &
                     + c * n(i,k,f,e,2) * (up(i,k,f,e,2) - u(i,j,k,e,2))  &
                     + c * n(i,k,f,e,3) * (up(i,k,f,e,3) - u(i,j,k,e,3))
        end do
        end do
      end do

      ! faces 5 and 6
      do f = 5, 6
        k = 1 + (np - 1) * (f - 5)
        do j = 1, np
        do i = 1, np
          c = a(i,j,f,e) / (2 * Ms(k) * Jd(i,j,k,e))
          v(i,j,k,e) = v(i,j,k,e)                                         &
                     + c * n(i,j,f,e,1) * (up(i,j,f,e,1) - u(i,j,k,e,1))  &
                     + c * n(i,j,f,e,2) * (up(i,j,f,e,2) - u(i,j,k,e,2))  &
                     + c * n(i,j,f,e,3) * (up(i,j,f,e,3) - u(i,j,k,e,3))
        end do
        end do
      end do

    end if

  end do

end subroutine TPO_Div_D_Gen_RWP
