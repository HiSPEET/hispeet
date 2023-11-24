!> summary:  Incompressible Navier-Stokes DG-SEM diffusion term (DC)
!> author:   Joerg Stiller
!> date:     2022/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetDiffusionTerm_C

  use C_Binding
  use TPO__Stress__3D
  use Mesh_Element__3D
  use Element_Face_Transfer_Buffer__3D

#ifdef __LIBXSMM__
  use LIBXSMM, only: LIBXSMM_XMMCall_ABC,  &
                     LIBXSMM_XMMDispatch,  &
                     LIBXSMM_BLASINT_KIND, &
                     LIBXSMM_DATATYPE_F32, &
                     LIBXSMM_DATATYPE_F64
#endif

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity on irregular (deformed) mesh
  !>
  !> The diffusion term `F_d` comprises the viscous contributions to the weak
  !> form of the momentum equation including the boundary conditions.
  !> To obtain the node values, apply the inverse mass matrix to `F_d`.

  module subroutine GetDiffusionTerm_C(this, v, vp, sp, F_d, bv_u, xout, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity (np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: vp(:,:,:,:,:)
    !< outer velocity traces (np,np,6,ne,3): vp = v⁺

    real(RNP), contiguous, intent(out) :: sp(:,:,:,:,:)
    !< outer viscous flux traces (np,np,6,ne,3): sp = s⁺ = -n⋅τ⁺,

    real(RNP), contiguous, intent(out) :: F_d(:,:,:,:,:)
    !< diffusion term (np,np,np,ne,3)

    class(BoundaryVariable_3D), optional, intent(in) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  vᵇ    in components 1:3
    !!   - Γᴼ :  τ_nn  in component    4, including penalty term
    !!
    !! ignored at extrapolated boundaries

    logical, optional, intent(in) :: xout
    !< ignore boundary values and extrapolate `∇v` [F]

    integer, optional, intent(in) :: form
    !< switch form of the approximated diffusion term `∇⋅τ`
    !!   - `0` :  `∇⋅ν[∇v + (∇v)ᵀ] + ∇[(μ-2/3ν)∇⋅v]`               (default)
    !!   - `1` :  `∇⋅ν[∇v + (∇v)ᵀ] - ∇[      ν ∇⋅v] =  ν∇²v`       (diffusion)
    !!   - `2` :  `∇⋅ν[∇v + (∇v)ᵀ] - ∇[     2ν ∇⋅v] = -ν∇×(∇×v)`   (rotation)

    ! local variables ..........................................................

    type(BoundaryVariable_3D), allocatable, save :: bv_s(:)
    ! viscous boundary flux

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr
    ! trace transfer buffers operators

    real(RNP), allocatable, save :: tr(:,:,:,:,:)
    ! traces of v and n⋅τ

    real(RNP) :: mu, nu
    integer   :: n_bound, n_elem, n_ghost, po
    integer   :: b

    ! initialization ...........................................................

    nu = this % nu_0

    if (present(form)) then
      select case(form)
      case(1)
        mu = -THIRD * nu
      case(2)
        mu = -THIRD * 4 * nu
    case default
        mu = this % mu_0
      end select
    else
      mu = this % mu_0
    end if

    n_bound = this % mesh % n_bound
    n_elem  = this % mesh % n_elem
    n_ghost = this % mesh % n_ghost
    po      = this % eop_v % po

    !$omp master
    allocate(bv_s(n_bound), tr(0:po,0:po,6,n_elem + n_ghost,6))
    buf_tr = ElementFaceTransferBuffer_3D(this%mesh, tr)
    do b = 1, n_bound
      if (scan(this % bc_v(b),'O') > 0) then
        call bv_s(b) % Create(this%mesh%boundary(b), po, nc = 3)
      end if
    end do
    !$omp end master
    !$omp barrier

    associate( mesh   => this % sem_v % mesh         &
             , Jd     => this % sem_v % metrics % Jd &
             , Ji     => this % sem_v % metrics % Ji &
             , a      => this % sem_v % metrics % a  &
             , n      => this % sem_v % metrics % n  &
             , eop    => this % eop_v                &
             , vm     => tr(:,:,:,:,1:3)             &
             , sm     => tr(:,:,:,:,4:6)             )

      ! element volume contributions ...........................................

      ! element integrals and local contributions to inner traces v⁻ and s⁻
      call TPO_Stress(eop%w, eop%D, Jd, Ji, n, mu, nu, v, F_d, vm, sm)

      ! transfer traces ........................................................

      call buf_tr % Transfer(mesh, tr, tag=1000) ! masters send inner traces
      call buf_tr % Merge(tr)                    ! ghosts  recv inner traces

      ! outer traces with BC ...................................................

      ! compute viscous boundary flux, where required
      do b = 1, mesh % n_bound
        if (scan(this % bc_v(b),'O') > 0) then
          associate(sb => bv_s(b) % val)
            call this % GetViscousBoundaryStress(b, v, sb, bv_u, xout, form)
          end associate
        end if
      end do

      call ConvertInnerToOuterTraces_3D(mesh, vm, vp)
      call this % ApplyEssentialBC(bv_u, vm, vp)       ! v⁺ = v(v⁻,vᵇ)  ∂Ω

      call ConvertInnerToOuterTraces_3D(mesh, sm, sp)
      call this % ApplyNaturalBC(bv_s, sm, sp)         ! s⁺ = s(s⁻,sᵇ)  ∂Ω

      ! fluxes .................................................................

      call AddViscousFluxes(mesh, eop, a, n, Ji, mu, nu, vm, sm, vp, sp, F_d )

      ! cleanup ................................................................

      !$omp master
      deallocate(bv_s, tr, buf_tr)
      !$omp end master

    end associate

  end subroutine GetDiffusionTerm_C

  !-----------------------------------------------------------------------------
  !> Compute & add viscous fluxes

  subroutine AddViscousFluxes(mesh, eop, a, n, Ji, mu, nu, vm, sm, vp, sp, F_d)

    class(Mesh_3D), intent(in) :: mesh
    !< mesh partition

    class(DG_ElementOperators_1D), intent(in) :: eop
    !< ID-DG element operators

    real(RNP), contiguous, intent(in) :: a(0:,0:,:,:)
    !< area coefficient at element faces (0:po,0:po,6,ne)

    real(RNP), contiguous, intent(in) :: n(0:,0:,:,:,:)
    !< unit normal vector at element faces (0:po,0:po,6,ne,3)

    real(RNP), contiguous, intent(in) :: Ji(0:,0:,0:,:,:,:)
    !< inverse Jacobian of element mapping (0:po,0:po,0:po,ne,3,3)

    real(RNP), intent(in) :: mu
    !< kinematic bulk viscosity

    real(RNP), intent(in) :: nu
    !< kinematic shear viscosity

    real(RNP), contiguous, intent(in) :: vm(0:,0:,:,:,:)
    !< v⁻ at element faces, including ghosts (0:po,0:po,6,ne+ng,3)

    real(RNP), contiguous, intent(in) :: sm(0:,0:,:,:,:)
    !< s⁻ = n⋅τ⁻ at element faces, including ghosts (0:po,0:po,6,ne+ng,3)

    real(RNP), contiguous, intent(in) :: vp(0:,0:,:,:,:)
    !< outer traces of velocity, v⁺ (0:po,0:po,6,ne,3)

    real(RNP), contiguous, intent(in) :: sp(0:,0:,:,:,:)
    !< outer traces of viscous stress vector, s⁺ = n⁺⋅τ⁺ (0:po,0:po,6,ne,3)

    real(RNP), contiguous, intent(inout) :: F_d(0:,0:,0:,:,:)
    !< diffusion term (0:po,0:po,0:po,ne,3)

    ! local variables ..........................................................

    type(C_FunPtr), save :: xmm

    ! target attribute required in call to LIBXSMM_XMMCall_ABC
    real(RNP), target, dimension(0:eop%po, 0:eop%po)       :: M_f, Ds_t
    real(RNP), target, dimension(0:eop%po, 0:eop%po, 3)    :: C0
    real(RNP), target, dimension(0:eop%po, 0:eop%po, 2, 3) :: Ct, Ds_Ct
    real(RNP), target, dimension(0:eop%po, 0:eop%po, 6, 3) :: Cn
    real(RNP) :: chi, cip, dv1, dv2, dv3, n1, n2, n3, tau(3,3)
    real(RNP) :: tmp
    integer   :: c, e, f, i, j, k

#ifdef __LIBXSMM__

    !$omp single
    select case(RNP)
    case(RSP)
      call DispatchKernel_RSP(xmm, eop%po)
    case(RDP)
      call DispatchKernel_RDP(xmm, eop%po)
    case default
      call Error( 'AddFluxes'                               &
                , 'Precision not supported'                 &
                , 'INS__Operator__3D:MP_GetDiffusionTerm_C' )
    end select
    !$omp end single

#endif

    associate(po => eop % po, Ms => eop % w, Ds => eop % D)

      ! auxiliaries ............................................................

      chi = mu - 2 * THIRD * nu

      ! face standard mass matrix
      do j = 0, po
      do i = 0, po
        M_f(i,j) = Ms(i) * Ms(j)
      end do
      end do

      Ds_t = transpose(Ds)

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east) ............................................

          do f = 1, 2

            i = (f - 1) * po

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do k = 0, po
            do j = 0, po

              n1 = n(j,k,f,e,1)
              n2 = n(j,k,f,e,2)
              n3 = n(j,k,f,e,3)

              dv1 = vm(j,k,f,e,1) - vp(j,k,f,e,1)
              dv2 = vm(j,k,f,e,2) - vp(j,k,f,e,2)
              dv3 = vm(j,k,f,e,3) - vp(j,k,f,e,3)

              tmp = chi * (n1 * dv1 + n2 * dv2 + n3 * dv3 )

              tau(1,1) = nu * 2 * n1 * dv1 + tmp
              tau(2,2) = nu * 2 * n2 * dv2 + tmp
              tau(3,3) = nu * 2 * n3 * dv3 + tmp

              tau(1,2) = nu * (n1 * dv2 + n2 * dv1)
              tau(2,3) = nu * (n2 * dv3 + n3 * dv2)
              tau(3,1) = nu * (n3 * dv1 + n1 * dv3)

              tau(2,1) = tau(1,2)
              tau(1,3) = tau(3,1)
              tau(3,2) = tau(2,3)

              tmp = HALF * M_f(j,k) * a(j,k,f,e)

              ! ξ-direction is "normal"
              do c = 1, 3

                Cn(j,k,f,c) = tmp * ( Ji(i,j,k,e,1,1) * tau(1,c) &
                                    + Ji(i,j,k,e,1,2) * tau(2,c) &
                                    + Ji(i,j,k,e,1,3) * tau(3,c) )

                Ct(j,k,1,c) = tmp * ( Ji(i,j,k,e,2,1) * tau(1,c) &
                                    + Ji(i,j,k,e,2,2) * tau(2,c) &
                                    + Ji(i,j,k,e,2,3) * tau(3,c) )

                Ct(k,j,2,c) = tmp * ( Ji(i,j,k,e,3,1) * tau(1,c) & ! first two
                                    + Ji(i,j,k,e,3,2) * tau(2,c) & ! directions
                                    + Ji(i,j,k,e,3,3) * tau(3,c) ) ! transposed

                C0(j,k,c)   = tmp * ( sm(j,k,f,e,c)               &
                                    - sp(j,k,f,e,c)               &
                                    - 2 * cip * ( n1 * tau(1,c)   &
                                                + n2 * tau(2,c)   &
                                                + n3 * tau(3,c) ) )

              end do

            end do
            end do

            select case(po)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))
              ! passing C addresses required by C interface
#else
              call TangentialDerivatives(po, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do k = 0, po
              do j = 0, po

                F_d(i,j,k,e,c) = F_d(i,j,k,e,c)  &
                               + Ds_Ct(j,k,1,c)  &
                               + Ds_Ct(k,j,2,c)  &
                               + C0(j,k,c)
              end do
              end do
            end do

          end do ! faces 1+2

          ! faces 3+4 (south + north) ..........................................

          do f = 3, 4

            j = (f - 3) * po

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do k = 0, po
            do i = 0, po

              n1 = n(i,k,f,e,1)
              n2 = n(i,k,f,e,2)
              n3 = n(i,k,f,e,3)

              dv1 = vm(i,k,f,e,1) - vp(i,k,f,e,1)
              dv2 = vm(i,k,f,e,2) - vp(i,k,f,e,2)
              dv3 = vm(i,k,f,e,3) - vp(i,k,f,e,3)

              tmp = chi * (n1 * dv1 + n2 * dv2 + n3 * dv3 )

              tau(1,1) = nu * 2 * n1 * dv1 + tmp
              tau(2,2) = nu * 2 * n2 * dv2 + tmp
              tau(3,3) = nu * 2 * n3 * dv3 + tmp

              tau(1,2) = nu * (n1 * dv2 + n2 * dv1)
              tau(2,3) = nu * (n2 * dv3 + n3 * dv2)
              tau(3,1) = nu * (n3 * dv1 + n1 * dv3)

              tau(2,1) = tau(1,2)
              tau(1,3) = tau(3,1)
              tau(3,2) = tau(2,3)

              tmp = HALF * M_f(i,k) * a(i,k,f,e)

              ! η-direction is "normal"
              do c = 1, 3

                Cn(i,k,f,c) = tmp * ( Ji(i,j,k,e,2,1) * tau(1,c) &
                                    + Ji(i,j,k,e,2,2) * tau(2,c) &
                                    + Ji(i,j,k,e,2,3) * tau(3,c) )

                Ct(i,k,1,c) = tmp * ( Ji(i,j,k,e,1,1) * tau(1,c) &
                                    + Ji(i,j,k,e,1,2) * tau(2,c) &
                                    + Ji(i,j,k,e,1,3) * tau(3,c) )

                Ct(k,i,2,c) = tmp * ( Ji(i,j,k,e,3,1) * tau(1,c) & ! first two
                                    + Ji(i,j,k,e,3,2) * tau(2,c) & ! directions
                                    + Ji(i,j,k,e,3,3) * tau(3,c) ) ! transposed

                C0(i,k,c)   = tmp * ( sm(i,k,f,e,c)               &
                                    - sp(i,k,f,e,c)               &
                                    - 2 * cip * ( n1 * tau(1,c)   &
                                                + n2 * tau(2,c)   &
                                                + n3 * tau(3,c) ) )

              end do

            end do
            end do

            select case(po)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))
#else
              call TangentialDerivatives(po, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do k = 0, po
              do i = 0, po

                F_d(i,j,k,e,c) = F_d(i,j,k,e,c)  &
                               + Ds_Ct(i,k,1,c)  &
                               + Ds_Ct(k,i,2,c)  &
                               + C0(i,k,c)
              end do
              end do
            end do

          end do ! faces 3+4

          ! faces 5+6 (top + bottom) ...........................................

          do f = 5, 6

            k = (f - 5) * po

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do j = 0, po
            do i = 0, po

              n1 = n(i,j,f,e,1)
              n2 = n(i,j,f,e,2)
              n3 = n(i,j,f,e,3)

              dv1 = vm(i,j,f,e,1) - vp(i,j,f,e,1)
              dv2 = vm(i,j,f,e,2) - vp(i,j,f,e,2)
              dv3 = vm(i,j,f,e,3) - vp(i,j,f,e,3)

              tmp = chi * (n1 * dv1 + n2 * dv2 + n3 * dv3 )

              tau(1,1) = nu * 2 * n1 * dv1 + tmp
              tau(2,2) = nu * 2 * n2 * dv2 + tmp
              tau(3,3) = nu * 2 * n3 * dv3 + tmp

              tau(1,2) = nu * (n1 * dv2 + n2 * dv1)
              tau(2,3) = nu * (n2 * dv3 + n3 * dv2)
              tau(3,1) = nu * (n3 * dv1 + n1 * dv3)

              tau(2,1) = tau(1,2)
              tau(1,3) = tau(3,1)
              tau(3,2) = tau(2,3)

              tmp = HALF * M_f(i,j) * a(i,j,f,e)

              ! ζ-direction is "normal"
              do c = 1, 3

                Cn(i,j,f,c) = tmp * ( Ji(i,j,k,e,3,1) * tau(1,c) &
                                    + Ji(i,j,k,e,3,2) * tau(2,c) &
                                    + Ji(i,j,k,e,3,3) * tau(3,c) )

                Ct(i,j,1,c) = tmp * ( Ji(i,j,k,e,1,1) * tau(1,c) &
                                    + Ji(i,j,k,e,1,2) * tau(2,c) &
                                    + Ji(i,j,k,e,1,3) * tau(3,c) )

                Ct(j,i,2,c) = tmp * ( Ji(i,j,k,e,2,1) * tau(1,c) & ! first two
                                    + Ji(i,j,k,e,2,2) * tau(2,c) & ! directions
                                    + Ji(i,j,k,e,2,3) * tau(3,c) ) ! transposed

                C0(i,j,c)   = tmp * ( sm(i,j,f,e,c)               &
                                    - sp(i,j,f,e,c)               &
                                    - 2 * cip * ( n1 * tau(1,c)   &
                                                + n2 * tau(2,c)   &
                                                + n3 * tau(3,c) ) )

              end do

            end do
            end do

            select case(po)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall_ABC(xmm, c_loc(Ds_t), c_loc(Ct), c_loc(Ds_Ct))
#else
              call TangentialDerivatives(po, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do j = 0, po
              do i = 0, po

                F_d(i,j,k,e,c) = F_d(i,j,k,e,c)  &
                               + Ds_Ct(i,j,1,c)  &
                               + Ds_Ct(j,i,2,c)  &
                               + C0(i,j,c)
              end do
              end do
            end do

          end do ! faces 5+6

         ! contributions to all points ...........................................

          do c = 1, 3
            do k = 0, po
            do j = 0, po
            do i = 0, po
              F_d(i,j,k,e,c) = F_d(i,j,k,e,c)           &
                             - ( Ds( 0,i) * Cn(j,k,1,c)  &
                               + Ds(po,i) * Cn(j,k,2,c)  &
                               + Ds( 0,j) * Cn(i,k,3,c)  &
                               + Ds(po,j) * Cn(i,k,4,c)  &
                               + Ds( 0,k) * Cn(i,j,5,c)  &
                               + Ds(po,k) * Cn(i,j,6,c)  &
                               )
            end do
            end do
            end do
          end do

        end associate
      end do
    end associate

  end subroutine AddViscousFluxes

  !-----------------------------------------------------------------------------
  !> Computation of tangential derivatives

  pure subroutine TangentialDerivatives(po, D, C, DC)
    integer,   intent(in)  :: po                 !< polynomial order
    real(RNP), intent(in)  :: D  (0:po,0:po)     !< diff matrix
    real(RNP), intent(in)  :: C  (0:po,0:po,2,3) !< tangential contributions
    real(RNP), intent(out) :: DC (0:po,0:po,2,3) !< derivatives of C

    integer :: i, j, l, m

    do l = 1, 3
      do j = 0, po
      do i = 0, po
        DC(i,j,1,l) = 0
        DC(i,j,2,l) = 0
        do m = 0, po
          DC(i,j,1,l) = DC(i,j,1,l) + D(m,i) * C(m,j,1,l)
          DC(i,j,2,l) = DC(i,j,2,l) + D(m,i) * C(m,j,2,l)
        end do
      end do
      end do
    end do

  end subroutine TangentialDerivatives

#ifdef __LIBXSMM__

  !-----------------------------------------------------------------------------
  !> Routine for dispatching a single precision XSMM kernel

  subroutine DispatchKernel_RSP(kernel, po, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: po
    integer(LIBXSMM_BLASINT_KIND), target, optional, intent(in) :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_FLOAT), target :: alpha = 1, beta = 0

    m = po + 1
    n = m * 6
    k = m

    call LIBXSMM_XMMDispatch( kernel, LIBXSMM_DATATYPE_F32       &
                            , m, n, k                            &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(alpha), c_loc(beta)          &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernel_RSP

  !-----------------------------------------------------------------------------
  !> Routine for dispatching a double precision XSMM kernel

  subroutine DispatchKernel_RDP(kernel, po, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: po
    integer(LIBXSMM_BLASINT_KIND), intent(in), optional, target :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_DOUBLE), target :: alpha = 1, beta = 0

    m = po + 1
    n = m * 6
    k = m

    call LIBXSMM_XMMDispatch( kernel, LIBXSMM_DATATYPE_F64       &
                            , m, n, k                            &
                            , c_loc(lda), c_loc(ldb), c_loc(ldc) &
                            , c_loc(alpha), c_loc(beta)          &
                            , c_loc(flags), c_loc(prefetch)      )

  end subroutine DispatchKernel_RDP

#endif

  !=============================================================================

end submodule MP_GetDiffusionTerm_C
