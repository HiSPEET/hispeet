!> summary:  3d curvilinear element stress operator kernel with LIBXSMM (DLCI)
!> author:   Joerg Stiller
!> date:     2022/02/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!>   - the grad-div-penalty version for incompressible Navier-Stokes
!>   - for consistency with surface terms, the penalty scales as χ = ζ - 2/3η
!===============================================================================

subroutine TPO_Stress_DLCI_XSMM_RWP(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)

  real(RWP), contiguous, intent(in) :: Ms(:)
  !< 1D diagonal standard mass matrix (np)

  real(RWP), contiguous, intent(in) :: Ds(:,:)
  !< 1D standard diff matrix (np,np)

  real(RWP), contiguous, intent(in) :: Jd(:,:,:,:)
  !< Jacobian determinant (np,np,np,ne)

  real(RWP), contiguous, intent(in) :: Ji(:,:,:,:,:,:)
  !< inverse Jacobian matrix  (np,np,np,ne,3,3)

  real(RWP), contiguous, intent(in) :: n(:,:,:,:,:)
  !< unit normal vector at element boundary faces (np,np,6,ne,3)

  real(RWP), intent(in) :: zeta
  !< bulk viscosity or modulus ζ

  real(RWP), intent(in) :: eta
  !< shear viscosity or modulus η

  real(RWP), contiguous, intent(in) :: v(:,:,:,:,:)
  !< velocity or displacement (np,np,np,ne,3)

  real(RWP), contiguous, intent(out) :: fv(:,:,:,:,:)
  !< element volume contribution to stress operator (np,np,np,ne,3)

  real(RWP), contiguous, optional, intent(out) :: vb(:,:,:,:,:)
  !< element boundary values of velocity or displacement (np,np,6,ne+ng,3),

  real(RWP), contiguous, optional, intent(out) :: sb(:,:,:,:,:)
  !< stress normal to element boundary fluxes, sb = n⁻⋅τ⁻ (np,np,6,ne+ng,3)

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: HALF      = 0.5_RWP
  real(RWP), parameter :: TWO_THIRD = 2.0_RWP / 3.0_RWP

  ! LIBXSMM kernel handles
  type(C_FunPtr), save :: xmm_grad(3), xmm_div(3)

  real(RWP), dimension(size(Ms), size(Ms))                 :: Ds_t
  real(RWP), dimension(size(Ms), size(Ms), size(Ms))       :: M, gdp
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 3)    :: w, z
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 3, 3) :: tau

  real(RWP) :: chi, dv(3)
  integer   :: c, e, f, i, j, k, np, ne
  logical   :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(v,4)

  chi = zeta - TWO_THIRD * eta

  get_traces = present(vb) .and. present(sb)

  Ds_t = transpose(Ds)

  ! XSMM kernel dispatch
  !$omp master
  select case(RWP)
  case(kind(1E0))
    call DispatchKernels_RSP(xmm_grad, xmm_div, np)
  case(kind(1D0))
    call DispatchKernels_RDP(xmm_grad, xmm_div, np)
  case default
    error stop '*** TPO_Stress_DLCI_XSMM_RWP: no match for RWP'
  end select
  !$omp end master
  !$omp barrier

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

    gdp = 0
    tau = 0

    do c = 1, 3

      ! standard derivatives of vc
      call LIBXSMM_XMMCall(xmm_grad(1), Ds, v(:,:,:,e,c), w(:,:,:,1))
      do k = 1, np
        call LIBXSMM_XMMCall(xmm_grad(2), v(:,:,k,e,c), Ds_t, w(:,:,k,2))
      end do
      call LIBXSMM_XMMCall(xmm_grad(3), v(:,:,:,e,c), Ds_t, w(:,:,:,3))

      ! tau = η(∇v + ∇vᵀ), gdp = χ∇⋅v
      do k = 1, np
      do j = 1, np
      do i = 1, np

        dv(1) = w(i,j,k,1) * Ji(i,j,k,e,1,1) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,1) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,1)

        dv(2) = w(i,j,k,1) * Ji(i,j,k,e,1,2) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,2) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,2)

        dv(3) = w(i,j,k,1) * Ji(i,j,k,e,1,3) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,3) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,3)

        tau(i,j,k,1,c) = tau(i,j,k,1,c) + eta * dv(1)
        tau(i,j,k,2,c) = tau(i,j,k,2,c) + eta * dv(2)
        tau(i,j,k,3,c) = tau(i,j,k,3,c) + eta * dv(3)

        tau(i,j,k,c,1) = tau(i,j,k,c,1) + eta * dv(1)
        tau(i,j,k,c,2) = tau(i,j,k,c,2) + eta * dv(2)
        tau(i,j,k,c,3) = tau(i,j,k,c,3) + eta * dv(3)

        gdp(i,j,k) = gdp(i,j,k) + chi * dv(c)

      end do
      end do
      end do

    end do

    ! volume integral ..........................................................

    ! w = M ΣᵣJ⁻¹(1,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * (tau(i,j,k,1,1) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,1,2) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,1,3) *  tau(i,j,k,3,1)               )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,1,1) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,1,2) * (tau(i,j,k,2,2) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,1,3) *  tau(i,j,k,3,2)               )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,1,1) *  tau(i,j,k,3,1)               &
                              + Ji(i,j,k,e,1,2) *  tau(i,j,k,3,2)               &
                              + Ji(i,j,k,e,1,3) * (tau(i,j,k,3,3) + gdp(i,j,k)) )
    end do
    end do
    end do

    ! z = [I x I x (Dˢ)ᵀ] w
    call LIBXSMM_XMMCall(xmm_div(1), Ds_t, w, z)

    ! w = M ΣᵣJ⁻¹(2,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * (tau(i,j,k,1,1) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,2,2) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,2,3) *  tau(i,j,k,3,1)               )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,2,1) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,2,2) * (tau(i,j,k,2,2) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,2,3) *  tau(i,j,k,3,2)               )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,2,1) *  tau(i,j,k,3,1)               &
                              + Ji(i,j,k,e,2,2) *  tau(i,j,k,3,2)               &
                              + Ji(i,j,k,e,2,3) * (tau(i,j,k,3,3) + gdp(i,j,k)) )
    end do
    end do
    end do

    ! z += [I x (Dˢ)ᵀ x I] w
    do c = 1, 3
      do k = 1, np
        call LIBXSMM_XMMCall(xmm_div(2), w(:,:,k,c), Ds, z(:,:,k,c))
      end do
    end do

    ! w = M ΣᵣJ⁻¹(3,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * (tau(i,j,k,1,1) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,3,2) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,3,3) *  tau(i,j,k,3,1)               )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,3,1) *  tau(i,j,k,2,1)               &
                              + Ji(i,j,k,e,3,2) * (tau(i,j,k,2,2) + gdp(i,j,k)) &
                              + Ji(i,j,k,e,3,3) *  tau(i,j,k,3,2)               )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,3,1) *  tau(i,j,k,3,1)               &
                              + Ji(i,j,k,e,3,2) *  tau(i,j,k,3,2)               &
                              + Ji(i,j,k,e,3,3) * (tau(i,j,k,3,3) + gdp(i,j,k)) )
    end do
    end do
    end do

    ! z += [(Dˢ)ᵀ x I x I] w
    do c = 1, 3
      call LIBXSMM_XMMCall(xmm_div(3), w(:,:,:,c), Ds, z(:,:,:,c))
    end do

    ! fv = -z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      fv(i,j,k,e,1) = -z(i,j,k,1)
      fv(i,j,k,e,2) = -z(i,j,k,2)
      fv(i,j,k,e,3) = -z(i,j,k,3)
    end do
    end do
    end do

    ! boundary values ..........................................................
    ! no penalty contribution here

    if (get_traces) then

      ! vb = v,  sb = n⋅τ  @ Γ₁ ∪ Γ₂, exploiting symmetry of τ

      i = 1
      do f = 1, 2
        do k = 1, np
        do j = 1, np

          vb(j,k,f,e,1) = v(i,j,k,e,1)
          vb(j,k,f,e,2) = v(i,j,k,e,2)
          vb(j,k,f,e,3) = v(i,j,k,e,3)

          sb(j,k,f,e,1) = n(j,k,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(j,k,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(j,k,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(j,k,f,e,2) = n(j,k,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(j,k,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(j,k,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(j,k,f,e,3) = n(j,k,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(j,k,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(j,k,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃

        end do
        end do
        i = np
      end do

      ! vb = v,  sb = n⋅τ  @ Γ₃ ∪ Γ₄

      j = 1
      do f = 3, 4
        do k = 1, np
        do i = 1, np

          vb(i,k,f,e,1) = v(i,j,k,e,1)
          vb(i,k,f,e,2) = v(i,j,k,e,2)
          vb(i,k,f,e,3) = v(i,j,k,e,3)

          sb(i,k,f,e,1) = n(i,k,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(i,k,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(i,k,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(i,k,f,e,2) = n(i,k,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(i,k,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(i,k,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(i,k,f,e,3) = n(i,k,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(i,k,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(i,k,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃
        end do
        end do
        j = np
      end do

      ! vb = v,  sb = n⋅τ  @ Γ₅ ∪ Γ₆

      k = 1
      do f = 5, 6
        do j = 1, np
        do i = 1, np

          vb(i,j,f,e,1) = v(i,j,k,e,1)
          vb(i,j,f,e,2) = v(i,j,k,e,2)
          vb(i,j,f,e,3) = v(i,j,k,e,3)

          sb(i,j,f,e,1) = n(i,j,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(i,j,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(i,j,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(i,j,f,e,2) = n(i,j,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(i,j,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(i,j,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(i,j,f,e,3) = n(i,j,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(i,j,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(i,j,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃

        end do
        end do
        k = np
      end do

    end if

  end do
  !$omp end do

contains

  !-----------------------------------------------------------------------------
  !> Routine for dispatching single precision XSMM kernels

  subroutine DispatchKernels_RSP( xmm_grad, xmm_div, np          &
                                , lda, ldb, ldc, flags, prefetch )

    type(C_FunPtr)               , intent(out) :: xmm_grad (3)
    type(C_FunPtr)               , intent(out) :: xmm_div  (3)
    integer(LIBXSMM_BLASINT_KIND), intent(in)  :: np
    integer(LIBXSMM_BLASINT_KIND), target, optional, intent(in) :: lda, ldb, ldc
    integer(C_INT)               , target, optional, intent(in) :: flags
    integer(C_INT)               , target, optional, intent(in) :: prefetch

    real(C_FLOAT), target :: ZERO, ONE
    integer(LIBXSMM_BLASINT_KIND) :: np2

    ZERO = 0
    ONE  = 1
    np2  = np * np

    call LIBXSMM_XMMDispatch( xmm_grad(1), LIBXSMM_DATATYPE_F32  &
                            , np, np2, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_grad(2), LIBXSMM_DATATYPE_F32  &
                            , np, np, np                         &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_grad(3), LIBXSMM_DATATYPE_F32  &
                            , np2, np, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(1), LIBXSMM_DATATYPE_F32   &
                            , np, np2*3, np                      &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(2), LIBXSMM_DATATYPE_F32   &
                            , np, np, np                         &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ONE)             &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(3), LIBXSMM_DATATYPE_F32   &
                            , np2, np, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ONE)             &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernels_RSP

  !-----------------------------------------------------------------------------
  !> Routine for dispatching double precision XSMM kernels

  subroutine DispatchKernels_RDP( xmm_grad, xmm_div, np          &
                                , lda, ldb, ldc, flags, prefetch )

    type(C_FunPtr)               , intent(out) :: xmm_grad (3)
    type(C_FunPtr)               , intent(out) :: xmm_div  (3)
    integer(LIBXSMM_BLASINT_KIND), intent(in)  :: np
    integer(LIBXSMM_BLASINT_KIND), target, optional, intent(in) :: lda, ldb, ldc
    integer(C_INT)               , target, optional, intent(in) :: flags
    integer(C_INT)               , target, optional, intent(in) :: prefetch

    real(C_DOUBLE), target :: ZERO, ONE
    integer(LIBXSMM_BLASINT_KIND) :: np2

    ZERO = 0
    ONE  = 1
    np2  = np * np

    call LIBXSMM_XMMDispatch( xmm_grad(1), LIBXSMM_DATATYPE_F64  &
                            , np, np2, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_grad(2), LIBXSMM_DATATYPE_F64  &
                            , np, np, np                         &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_grad(3), LIBXSMM_DATATYPE_F64  &
                            , np2, np, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(1), LIBXSMM_DATATYPE_F64   &
                            , np, np2*3, np                      &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ZERO)            &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(2), LIBXSMM_DATATYPE_F64   &
                            , np, np, np                         &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ONE)             &
                            , c_loc(flags), c_loc(prefetch)      )

    call LIBXSMM_XMMDispatch( xmm_div(3), LIBXSMM_DATATYPE_F64   &
                            , np2, np, np                        &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(ONE), c_loc(ONE)             &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernels_RDP

  !-----------------------------------------------------------------------------

end subroutine TPO_Stress_DLCI_XSMM_RWP
