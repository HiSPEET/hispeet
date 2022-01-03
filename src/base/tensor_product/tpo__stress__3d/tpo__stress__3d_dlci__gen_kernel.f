!> summary:  3d generic curvilinear element stress operator
!> author:   Joerg Stiller
!> date:     2021/12/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Stress_DLCI_Gen_RWP(Ms, Ds, Jd, Ji, n, eta, chi, v, &
                                   fv, v_b, tau_b                  )

  real(RWP), contiguous, intent(in) :: Ms(:)
  !< 1D diagonal standard mass matrix (np)

  real(RWP), contiguous, intent(in)  :: Ds(:,:)
  !< 1D standard diff matrix (np,np)

  real(RWP), contiguous, intent(in)  :: Jd(:,:,:,:)
  !< Jacobian determinant (np,np,np,ne)

  real(RWP), contiguous, intent(in)  :: Ji(:,:,:,:,:,:)
  !< inverse Jacobian matrix  (np,np,np,ne,3,3)

  real(RWP), contiguous, intent(in)  :: n(:,:,:,:,:)
  !< unit normal vector at element boundary faces (np,np,6,ne,3)

  real(RWP), intent(in)  :: eta
  !< shear viscosity or modulus η

  real(RWP), intent(in)  :: chi
  !< normalized first Lamé parameter, χ = λ/η; Stokes hypothesis: χ = -2/3

  real(RWP), contiguous, intent(in)  :: v(:,:,:,:,:)
  !< velocity or displacement (np,np,np,ne,3)

  real(RWP), contiguous, intent(out) :: fv(:,:,:,:,:)
  !< element volume contribution to stress operator (np,np,np,ne,3)

  real(RWP), contiguous, intent(out) :: v_b(:,:,:,:,:)
  !< element boundary values of velocity or displacement (np,np,6,ne,3)

  real(RWP), contiguous, intent(out) :: tau_b(:,:,:,:,:)
  !< stress normal to element boundary fluxes, tau_b = n⁻⋅τ⁻ (np,np,6,ne,3)

  !---------------------------------------------------------------------------
  ! local variables

  ! linear index of deformation tensor
  integer, parameter :: l(3,3) = reshape( [ 1, 2, 3,         &
                                            2, 4, 5,         &
                                            3, 5, 6 ], [3,3] )

  real(RWP), parameter :: HALF = 0.5_RWP

  real(RWP), dimension(size(Ms), size(Ms), size(Ms))    :: M
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 3) :: w, z
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 6) :: tau

  real(RWP) :: dv(3), def_v(6), chi_div_v
  real(RWP) :: tmp
  integer   :: np, ne
  integer   :: c, e, f, i, j, k, p
  integer   :: l1c, l2c, l3c, lc1, lc2, lc3

  !---------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(v,4)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! element mass matrix ......................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) = Ms(k) * Ms(j) * Ms(i) * Jd(i,j,k,e)
    end do
    end do
    end do

    ! viscous stress tensor ....................................................

    tau = 0

    do c = 1, 3

      ! w₁ = [I x I x Dˢ] v(:,:,:,e,c)
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(i,p) * v(p,j,k,e,c)
        end do
        w(i,j,k,1) = tmp
      end do
      end do
      end do

      ! w₂ = [I x Dˢ x I] v(:,:,:,e,c)
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(j,p) * v(i,p,k,e,c)
        end do
        w(i,j,k,2) = tmp
      end do
      end do
      end do

      ! w₃ = [Dˢ x I x I] v(:,:,:,e,c)
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(k,p) * v(i,j,p,e,c)
        end do
        w(i,j,k,3) = tmp
      end do
      end do
      end do

      ! tau = ∇v + ∇vᵀ

      l1c = l(1,c)
      l2c = l(2,c)
      l3c = l(3,c)
      lc1 = l(c,1)
      lc2 = l(c,2)
      lc3 = l(c,3)

      do k = 1, np
      do j = 1, np
      do i = 1, np

        dv(1) = Ji(i,j,k,e,1,1) * w(i,j,k,1) &
              + Ji(i,j,k,e,2,1) * w(i,j,k,2) &
              + Ji(i,j,k,e,3,1) * w(i,j,k,3)

        dv(2) = Ji(i,j,k,e,1,2) * w(i,j,k,1) &
              + Ji(i,j,k,e,2,2) * w(i,j,k,2) &
              + Ji(i,j,k,e,3,2) * w(i,j,k,3)

        dv(3) = Ji(i,j,k,e,1,3) * w(i,j,k,1) &
              + Ji(i,j,k,e,2,3) * w(i,j,k,2) &
              + Ji(i,j,k,e,3,3) * w(i,j,k,3)

        tau(i,j,k,l1c) = tau(i,j,k,l1c) + dv(1)
        tau(i,j,k,l2c) = tau(i,j,k,l2c) + dv(2)
        tau(i,j,k,l3c) = tau(i,j,k,l3c) + dv(3)

        tau(i,j,k,lc1) = tau(i,j,k,lc1) + dv(1)
        tau(i,j,k,lc2) = tau(i,j,k,lc2) + dv(2)
        tau(i,j,k,lc3) = tau(i,j,k,lc3) + dv(3)

      end do
      end do
      end do

    end do

    ! assemble stress tensor
    do k = 1, np
    do j = 1, np
    do i = 1, np

      def_v(1) = tau(i,j,k,1) !  ∂₁v₁ + ∂₁v₁
      def_v(2) = tau(i,j,k,2) !  ∂₁v₂ + ∂₂v₁
      def_v(3) = tau(i,j,k,3) !  ∂₁v₃ + ∂₃v₁
      def_v(4) = tau(i,j,k,4) !  ∂₂v₂ + ∂₂v₂
      def_v(5) = tau(i,j,k,5) !  ∂₂v₃ + ∂₃v₂
      def_v(6) = tau(i,j,k,6) !  ∂₃v₃ + ∂₃v₃

      chi_div_v = chi * HALF * (def_v(1) + def_v(4) + def_v(6)) !  χ ∇⋅v

      tau(i,j,k,1) = eta * (def_v(1) + chi_div_v) !  τ₁₁
      tau(i,j,k,2) = eta *  def_v(2)              !  τ₁₂ = τ₂₁
      tau(i,j,k,3) = eta *  def_v(3)              !  τ₁₃ = τ₃₁
      tau(i,j,k,4) = eta * (def_v(4) + chi_div_v) !  τ₂₂
      tau(i,j,k,5) = eta *  def_v(5)              !  τ₂₃ = τ₃₂
      tau(i,j,k,6) = eta * (def_v(6) + chi_div_v) !  τ₃₃

    end do
    end do
    end do

    ! volume integral ..........................................................

    ! w = M ΣᵣJ⁻¹(1,r)⋅τ(r,1:3)
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,1) & ! τ₁₁
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,2) & ! τ₂₁
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,3) ) ! τ₃₁

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,2) & ! τ₁₂
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,4) & ! τ₂₂
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,5) ) ! τ₃₂

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,3) & ! τ₁₃
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,5) & ! τ₂₃
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,6) ) ! τ₃₃
    end do
    end do
    end do

    ! z = [I x I x (Dˢ)ᵀ] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k,1) = 0
      z(i,j,k,2) = 0
      z(i,j,k,3) = 0
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,i) * w(p,j,k,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,i) * w(p,j,k,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,i) * w(p,j,k,3)
      end do
    end do
    end do
    end do

    ! w = M ΣᵣJ⁻¹(2,r)⋅τ(r,1:3)
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,1) & ! τ₁₁
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,2) & ! τ₂₁
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,3) ) ! τ₃₁

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,2) & ! τ₁₂
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,4) & ! τ₂₂
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,5) ) ! τ₃₂

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,3) & ! τ₁₃
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,5) & ! τ₂₃
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,6) ) ! τ₃₃
    end do
    end do
    end do

    ! z += [I x (Dˢ)ᵀ x I] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,j) * w(i,p,k,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,j) * w(i,p,k,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,j) * w(i,p,k,3)
      end do
    end do
    end do
    end do

    ! w = M ΣᵣJ⁻¹(3,r)⋅τ(r,1:3)
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,1) & ! τ₁₁
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,2) & ! τ₂₁
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,3) ) ! τ₃₁

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,2) & ! τ₁₂
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,4) & ! τ₂₂
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,5) ) ! τ₃₂

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,3) & ! τ₁₃
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,5) & ! τ₂₃
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,6) ) ! τ₃₃
    end do
    end do
    end do

    ! z += [(Dˢ)ᵀ x I x I] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,k) * w(i,j,p,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,k) * w(i,j,p,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,k) * w(i,j,p,3)
      end do
    end do
    end do
    end do

    ! fv = z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      fv(i,j,k,e,1) = z(i,j,k,1)
      fv(i,j,k,e,2) = z(i,j,k,2)
      fv(i,j,k,e,3) = z(i,j,k,3)
    end do
    end do
    end do

    ! boundary values ..........................................................

    ! v_b = v,  tau_b = n⋅τ  @ Γ₁ ∪ Γ₂

    i = 1
    do f = 1, 2
      do k = 1, np
      do j = 1, np

        v_b(j,k,f,e,1) = v(i,j,k,e,1)
        v_b(j,k,f,e,2) = v(i,j,k,e,2)
        v_b(j,k,f,e,3) = v(i,j,k,e,3)

        tau_b(j,k,f,e,1) = n(j,k,f,e,1) * tau(i,j,k,1) & ! n₁ τ₁₁
                         + n(j,k,f,e,2) * tau(i,j,k,2) & ! n₂ τ₂₁
                         + n(j,k,f,e,3) * tau(i,j,k,3)   ! n₃ τ₃₁

        tau_b(j,k,f,e,2) = n(j,k,f,e,1) * tau(i,j,k,2) & ! n₁ τ₁₂
                         + n(j,k,f,e,2) * tau(i,j,k,4) & ! n₂ τ₂₂
                         + n(j,k,f,e,3) * tau(i,j,k,5)   ! n₃ τ₃₂

        tau_b(j,k,f,e,3) = n(j,k,f,e,1) * tau(i,j,k,3) & ! n₁ τ₁₃
                         + n(j,k,f,e,2) * tau(i,j,k,5) & ! n₂ τ₂₃
                         + n(j,k,f,e,3) * tau(i,j,k,6)   ! n₃ τ₃₃

      end do
      end do
      i = np
    end do

    ! v_b = v,  tau_b = n⋅τ  @ Γ₃ ∪ Γ₄

    j = 1
    do f = 3, 4
      do k = 1, np
      do i = 1, np

        v_b(i,k,f,e,1) = v(i,j,k,e,1)
        v_b(i,k,f,e,2) = v(i,j,k,e,2)
        v_b(i,k,f,e,3) = v(i,j,k,e,3)

        tau_b(i,k,f,e,1) = n(i,k,f,e,1) * tau(i,j,k,1) & ! n₁ τ₁₁
                         + n(i,k,f,e,2) * tau(i,j,k,2) & ! n₂ τ₂₁
                         + n(i,k,f,e,3) * tau(i,j,k,3)   ! n₃ τ₃₁

        tau_b(i,k,f,e,2) = n(i,k,f,e,1) * tau(i,j,k,2) & ! n₁ τ₁₂
                         + n(i,k,f,e,2) * tau(i,j,k,4) & ! n₂ τ₂₂
                         + n(i,k,f,e,3) * tau(i,j,k,5)   ! n₃ τ₃₂

        tau_b(i,k,f,e,3) = n(i,k,f,e,1) * tau(i,j,k,3) & ! n₁ τ₁₃
                         + n(i,k,f,e,2) * tau(i,j,k,5) & ! n₂ τ₂₃
                         + n(i,k,f,e,3) * tau(i,j,k,6)   ! n₃ τ₃₃
      end do
      end do
      j = np
    end do

    ! v_b = v,  tau_b = n⋅τ  @ Γ₅ ∪ Γ₆

    k = 1
    do f = 5, 6
      do j = 1, np
      do i = 1, np

        v_b(i,j,f,e,1) = v(i,j,k,e,1)
        v_b(i,j,f,e,2) = v(i,j,k,e,2)
        v_b(i,j,f,e,3) = v(i,j,k,e,3)

        tau_b(i,j,f,e,1) = n(i,j,f,e,1) * tau(i,j,k,1) & ! n₁ τ₁₁
                         + n(i,j,f,e,2) * tau(i,j,k,2) & ! n₂ τ₂₁
                         + n(i,j,f,e,3) * tau(i,j,k,3)   ! n₃ τ₃₁

        tau_b(i,j,f,e,2) = n(i,j,f,e,1) * tau(i,j,k,2) & ! n₁ τ₁₂
                         + n(i,j,f,e,2) * tau(i,j,k,4) & ! n₂ τ₂₂
                         + n(i,j,f,e,3) * tau(i,j,k,5)   ! n₃ τ₃₂

        tau_b(i,j,f,e,3) = n(i,j,f,e,1) * tau(i,j,k,3) & ! n₁ τ₁₃
                         + n(i,j,f,e,2) * tau(i,j,k,5) & ! n₂ τ₂₃
                         + n(i,j,f,e,3) * tau(i,j,k,6)   ! n₃ τ₃₃

      end do
      end do
      k = np
    end do

  end do
  !$omp end do

end subroutine TPO_Stress_DLCI_Gen_RWP
