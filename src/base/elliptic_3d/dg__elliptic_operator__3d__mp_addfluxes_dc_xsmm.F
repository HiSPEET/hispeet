!> summary:  3D DG elliptic operator: application of fluxes with deformed mesh
!>           and constant diffusivity using LIBXSMM
!> author:   Joerg Stiller
!> date:     2021/11/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D:MP_Eval_DC) MP_AddFluxes_DC_XSMM

#ifdef __LIBXSMM__

  use, intrinsic ::  ISO_C_Binding

  use LIBXSMM, only: LIBXSMM_XMMCall_ABC,  &
                     LIBXSMM_XMMDispatch,  &
                     LIBXSMM_BLASINT_KIND, &
                     LIBXSMM_DATATYPE_F32, &
                     LIBXSMM_DATATYPE_F64

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes: LIBXSMM version

  module subroutine AddFluxes_DC_XSMM(this, hom_bc, a, Ji_n, nu, tr, r, f)

    ! arguments ................................................................

    class(DG_EllipticOperator_3D), intent(in) :: this

    logical,   intent(in)    :: hom_bc            !< F/T for in/homogeneous BC
    real(RNP), intent(in)    :: a(0:,0:,:,:)      !< area coeff @ element faces
    real(RNP), intent(in)    :: Ji_n(0:,0:,:,:,:) !< J⁻¹⋅n      @ element faces
    real(RNP), intent(in)    :: nu                !< diffusivity, ν = νᵖ
    real(RNP), intent(in)    :: tr(0:,0:,:,:,:)   !< traces of u, q_n
    real(RNP), intent(inout) :: r(0:,0:,0:,:)     !< result
    real(RNP), intent(in), optional :: f(0:,0:,0:,:) !< RHS

    contiguous :: a, Ji_n, tr, r, f

    ! local variables ..........................................................

    ! LIBXSMM kernel handle
    type(C_FunPtr), save :: xmm

    ! target attribute required in call to LIBXSMM_XMMCall_ABC
    real(RNP), target, dimension(0:this%eop%po, 0:this%eop%po) :: Mf, Ds_t
    real(RNP), target, dimension(0:this%eop%po, 0:this%eop%po) :: jmp_u, avg_q
    real(RNP), target, dimension(0:this%eop%po, 0:this%eop%po, 2) :: Ct, Ds_Ct
    real(RNP), target, dimension(0:this%eop%po, 0:this%eop%po, 6) :: Cn
    real(RNP) :: mu_nu, tmp
    integer   :: e, i, j, k, m
    logical   :: residual, struct

    associate( eop  => this % eop        &
             , P    => this % eop % po   &
             , Ms   => this % eop % w    &
             , mesh => this % sem % mesh )

      ! auxiliaries ............................................................

      !$omp single

      ! dispatch XSMM kernel
      select case(RNP)
      case(RSP)
        call DispatchKernel_RSP(xmm, P)
      case(RDP)
        call DispatchKernel_RDP(xmm, P)
      case default
        call Error( 'AddFluxes_XSMM'                                  &
                  , 'Precision not supported'                         &
                  , 'DG__Elliptic_Operator__3D:MP_AddFLuxes_DC_XSMM' )
      end select
      !$omp end single

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      Ds_t = transpose(eop % D)

      residual = present(f)
      struct   = mesh % structured

      ! add fluxes .............................................................

      !$omp do private(e)
      do e = 1, mesh % n_elem_active
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east)  . . . . . . . . . . . . . . . . . . . . . .

          do m = 1, 2

            ! jmp_u = n·[u], avg_q = n·{ν∇u}
            call this % GetElementBoundaryFluxes( element, struct, hom_bc &
                                                , e, m, tr, jmp_u, avg_q  )

            mu_nu = eop % PenaltyFactor(element % geometry % dx_m(m)) * nu
            i = (m-1) * P

            do k = 0, P
            do j = 0, P
              tmp = Mf(j,k) * HALF * nu * a(j,k,m,e) * jmp_u(j,k)
              Cn(j,k,m) = tmp * Ji_n(j,k,m,e,1)
              Ct(j,k,1) = tmp * Ji_n(j,k,m,e,2) ! Ct(:,:,2) is transposed for
              Ct(k,j,2) = tmp * Ji_n(j,k,m,e,3) ! fused evaluation of derivative
            end do
            end do

            ! derivatives of tangential jump contributions
            call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))
            ! passing C addresses required by C interface

            ! add flux contributions, accounting for transposition of Ct(:,:,2)
            do k = 0, P
            do j = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(j,k) * a(j,k,m,e)               &
                           * (avg_q(j,k) - mu_nu * jmp_u(j,k))  &
                           + Ds_Ct(j,k,1)                       &
                           + Ds_Ct(k,j,2)                       &
                           )
            end do
            end do

          end do

          ! faces 3+4 (south + north)  . . . . . . . . . . . . . . . . . . . . .

          do m = 3, 4

            call this % GetElementBoundaryFluxes( element, struct, hom_bc &
                                                , e, m, tr, jmp_u, avg_q  )

            mu_nu = eop % PenaltyFactor(element % geometry % dx_m(m)) * nu
            j = (m-3) * P

            do k = 0, P
            do i = 0, P
              tmp = Mf(i,k) * HALF * nu * a(i,k,m,e) * jmp_u(i,k)
              Cn(i,k,m) = tmp * Ji_n(i,k,m,e,2)
              Ct(i,k,1) = tmp * Ji_n(i,k,m,e,1)
              Ct(k,i,2) = tmp * Ji_n(i,k,m,e,3)
            end do
            end do

            call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))

            do k = 0, P
            do i = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(i,k) * a(i,k,m,e)               &
                           * (avg_q(i,k) - mu_nu * jmp_u(i,k))  &
                           + Ds_Ct(i,k,1)                       &
                           + Ds_Ct(k,i,2)                       &
                           )
            end do
            end do

          end do

          ! faces 5+6 (bottom + top) . . . . . . . . . . . . . . . . . . . . . .

          do m = 5, 6

            call this % GetElementBoundaryFluxes( element, struct, hom_bc &
                                                , e, m, tr, jmp_u, avg_q  )

            mu_nu = eop % PenaltyFactor(element % geometry % dx_m(m)) * nu
            k = (m-5) * P

            do j = 0, P
            do i = 0, P
              tmp = Mf(i,j) * HALF * nu * a(i,j,m,e) * jmp_u(i,j)
              Cn(i,j,m) = tmp * Ji_n(i,j,m,e,3)
              Ct(i,j,1) = tmp * Ji_n(i,j,m,e,1)
              Ct(j,i,2) = tmp * Ji_n(i,j,m,e,2)
            end do
            end do

            call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))

            do j = 0, P
            do i = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(i,j) * a(i,j,m,e)               &
                           * (avg_q(i,j) - mu_nu * jmp_u(i,j))  &
                           + Ds_Ct(i,j,1)                       &
                           + Ds_Ct(j,i,2)                       &
                           )
            end do
            end do

          end do

          if (residual) then
            ! r = f - Au
            do k = 0, P
            do j = 0, P
            do i = 0, P
              r(i,j,k,e) = Ds_t(i,0) * Cn(j,k,1)  &
                         + Ds_t(i,P) * Cn(j,k,2)  &
                         + Ds_t(j,0) * Cn(i,k,3)  &
                         + Ds_t(j,P) * Cn(i,k,4)  &
                         + Ds_t(k,0) * Cn(i,j,5)  &
                         + Ds_t(k,P) * Cn(i,j,6)  &
                         + f(i,j,k,e)             &
                         - r(i,j,k,e)
            end do
            end do
            end do
          else
            ! r = Au
            do k = 0, P
            do j = 0, P
            do i = 0, P
              r(i,j,k,e) = r(i,j,k,e)               &
                         - ( Ds_t(i,0) * Cn(j,k,1)  &
                           + Ds_t(i,P) * Cn(j,k,2)  &
                           + Ds_t(j,0) * Cn(i,k,3)  &
                           + Ds_t(j,P) * Cn(i,k,4)  &
                           + Ds_t(k,0) * Cn(i,j,5)  &
                           + Ds_t(k,P) * Cn(i,j,6)  &
                           )
            end do
            end do
            end do
          end if

        end associate

      end do
    end associate

  end subroutine AddFluxes_DC_XSMM

  !-----------------------------------------------------------------------------
  !> Routine for dispatching a single precision XSMM kernel

  subroutine DispatchKernel_RSP(kernel, P, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: P
    integer(LIBXSMM_BLASINT_KIND), target, optional, intent(in) :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_FLOAT), target :: alpha = 1, beta = 0

    m = P + 1
    n = m * 2
    k = m

    call LIBXSMM_XMMDispatch( kernel, LIBXSMM_DATATYPE_F32       &
                            , m, n, k                            &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(alpha), c_loc(beta)          &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernel_RSP

  !-----------------------------------------------------------------------------
  !> Routine for dispatching a double precision XSMM kernel

  subroutine DispatchKernel_RDP(kernel, P, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: P
    integer(LIBXSMM_BLASINT_KIND), intent(in), optional, target :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_DOUBLE), target :: alpha = 1, beta = 0

    m = P + 1
    n = m * 2
    k = m

    call LIBXSMM_XMMDispatch( kernel, LIBXSMM_DATATYPE_F64       &
                            , m, n, k                            &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(alpha), c_loc(beta)          &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernel_RDP

  !=============================================================================

#endif

end submodule MP_AddFluxes_DC_XSMM
