!> summary:  Incompressible Navier-Stokes DG-SEM diffusion term (DC)
!> author:   Joerg Stiller
!> date:     2022/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetDiffusionTerm_DC

#ifdef __LIBXSMM__
  use LIBXSMM
#endif

  use C_Binding
  use TPO__Stress__3D
  use Mesh_Element__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity on irregular (deformed) mesh
  !>
  !> The diffusion term `F_d` comprises the viscous contributions to the weak
  !> form of the momentum equation including the boundary conditions.
  !> For homogeneous conditions pass `vp = sp = 0`.
  !> To obtain the node values, apply the inverse mass matrix to `F_d`.

  module subroutine GetDiffusionTerm_DC(this, v, vp, sp, F_d)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity (np,np,np,ne,3)

    real(RNP), contiguous, intent(inout) :: vp(:,:,:,:,:)
    !< exterior velocity traces (np,np,6,ne,3)
    !<   - in:  v   at Dirichlet faces, nn⋅v at free-slip faces, undefined else
    !<   - out: v⁺  at all element faces

    real(RNP), contiguous, intent(inout) :: sp(:,:,:,:,:)
    !< viscous flux traces (np,np,6,ne,3)
    !<   - in:  s  =  n⋅τ   at free-slip or traction faces, undefined else
    !<   - out: s⁺ = -n⋅τ⁺  at all element faces, since n⁺ = -n

    real(RNP), contiguous, intent(out) :: F_d(:,:,:,:,:)
    !< diffusion term (np,np,np,ne,3)

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of v and n⋅τ

    associate( mesh   => this % sem_v % mesh         &
             , Jd     => this % sem_v % metrics % Jd &
             , Ji     => this % sem_v % metrics % Ji &
             , a      => this % sem_v % metrics % a  &
             , n      => this % sem_v % metrics % n  &
             , eop    => this % eop_v                &
             , bc     => this % bc_v                 &
             , mu     => this % mu_0                 &
             , nu     => this % nu_0                 )

      ! initialization .........................................................

      ! workspace and operators
      !$omp master
      allocate(tr(0:eop%po, 0:eop%po, 6, mesh%n_elem + mesh%n_ghost,6))
      buf_tr = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! element volume contributions ...........................................

      call TPO_Stress( eop%w, eop%D, Jd, Ji, n, mu, nu, v &
                     , fv = F_d                           &
                     , vb = tr(:,:,:,:,1:3)               &
                     , sb = tr(:,:,:,:,4:6)               )

      ! transfer traces ........................................................

      call buf_tr % Transfer(mesh, tr, tag=1000)
      call buf_tr % Merge(tr)

      ! add fluxes .............................................................

      call AddViscousFluxes( eop%po, mesh, eop, bc &
                           , a, n, Ji, mu, nu      &
                           , vm = tr(:,:,:,:,1:3)  &
                           , sm = tr(:,:,:,:,4:6)  &
                           , vp = vp               &
                           , sp = sp               &
                           , F_d = F_d             )

      ! cleanup ................................................................

      !$omp master
      deallocate(tr, buf_tr)
      !$omp end master

    end associate

  end subroutine GetDiffusionTerm_DC

  !-----------------------------------------------------------------------------
  !> Compute & add viscous fluxes

  subroutine AddViscousFluxes( P, mesh, eop, bc, a, n, Ji  &
                             , mu, nu, vm, sm, vp, sp, F_d )

    integer, intent(in) :: P
    !< polynomial order

    class(Mesh_3D), intent(in) :: mesh
    !< mesh partition

    class(DG_ElementOperators_1D), intent(in) :: eop
    !< ID-DG element operators

    character, intent(in) :: bc(:)
    !< boundary conditions

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

    real(RNP), contiguous, intent(inout) :: vp(0:,0:,:,:,:)
    !< exterior velocity traces at element faces (0:po,0:po,6,ne,3)
    !<   - in:  v   at Dirichlet faces, nn⋅v at free-slip faces, undefined else
    !<   - out: v⁺  at all element faces

    real(RNP), contiguous, intent(inout) :: sp(0:,0:,:,:,:)
    !< vicous stress vector at element faces (0:po,0:po,6,ne,3)
    !<   - in:  s  = n⋅τ   at free-slip or traction faces, undefined else
    !<   - out: s⁺ = n⋅τ⁺  at all element faces

    real(RNP), contiguous, intent(inout) :: F_d(0:,0:,0:,:,:)
    !< diffusion term (0:po,0:po,0:po,ne,3)

    ! local variables ..........................................................

    type(C_FunPtr), save :: xmm

    real(RNP), dimension(0:P, 0:P)       :: M_f, Ds_t
    real(RNP), dimension(0:P, 0:P, 3)    :: C0
    real(RNP), dimension(0:P, 0:P, 2, 3) :: Ct, Ds_Ct
    real(RNP), dimension(0:P, 0:P, 6, 3) :: Cn
    real(RNP) :: chi, cip, dv1, dv2, dv3, n1, n2, n3, tau(3,3)
    real(RNP) :: tmp
    integer   :: c, e, f, i, j, k
    logical   :: struct

#ifdef __LIBXSMM__

    !$omp single
    select case(RNP)
    case(RSP)
      call DispatchKernel_RSP(xmm, P)
    case(RDP)
      call DispatchKernel_RDP(xmm, P)
    case default
      call Error( 'AddFluxes'                                &
                , 'Precision not supported'                  &
                , 'INS__Operator__3D:MP_GetDiffusionTerm_DC' )
    end select
    !$omp end single

#endif

    associate(Ms => eop % w, Ds => eop % D)

      ! auxiliaries ............................................................

      struct  = mesh % structured

      chi = mu - 2 * THIRD * nu

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        M_f(i,j) = Ms(i) * Ms(j)
      end do
      end do

      Ds_t = transpose(Ds)

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east) ............................................

          do f = 1, 2

            i = (f - 1) * P

            call GetExteriorTraces(P, element, struct, bc, e,f,n, vm,sm, vp,sp)

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do k = 0, P
            do j = 0, P

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

            select case(P)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall(xmm, Ds_t, Ct, Ds_Ct)
#else
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do k = 0, P
              do j = 0, P

                F_d(i,j,k,e,c) = F_d(i,j,k,e,c)  &
                               + Ds_Ct(j,k,1,c)  &
                               + Ds_Ct(k,j,2,c)  &
                               + C0(j,k,c)
              end do
              end do
            end do

          end do ! faces 1+2

          ! faces 3+4 (south + north) ............................................

          do f = 3, 4

            j = (f - 3) * P

            call GetExteriorTraces(P, element, struct, bc, e,f,n, vm,sm, vp,sp)

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do k = 0, P
            do i = 0, P

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

              do c = 1, 3

                Cn(i,k,f,c) = tmp * ( Ji(i,j,k,e,1,1) * tau(1,c) &
                                    + Ji(i,j,k,e,1,2) * tau(2,c) &
                                    + Ji(i,j,k,e,1,3) * tau(3,c) )

                Ct(i,k,1,c) = tmp * ( Ji(i,j,k,e,2,1) * tau(1,c) &
                                    + Ji(i,j,k,e,2,2) * tau(2,c) &
                                    + Ji(i,j,k,e,2,3) * tau(3,c) )

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

            select case(P)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall(xmm, Ds_t, Ct, Ds_Ct)
#else
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do k = 0, P
              do i = 0, P

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

            k = (f - 5) * P

            call GetExteriorTraces(P, element, struct, bc, e,f,n, vm,sm, vp,sp)

            cip = eop % PenaltyFactor(element % geometry % dx_m(f))

            do j = 0, P
            do i = 0, P

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

              do c = 1, 3

                Cn(i,j,f,c) = tmp * ( Ji(i,j,k,e,1,1) * tau(1,c) &
                                    + Ji(i,j,k,e,1,2) * tau(2,c) &
                                    + Ji(i,j,k,e,1,3) * tau(3,c) )

                Ct(i,j,1,c) = tmp * ( Ji(i,j,k,e,2,1) * tau(1,c) &
                                    + Ji(i,j,k,e,2,2) * tau(2,c) &
                                    + Ji(i,j,k,e,2,3) * tau(3,c) )

                Ct(j,i,2,c) = tmp * ( Ji(i,j,k,e,3,1) * tau(1,c) & ! first two
                                    + Ji(i,j,k,e,3,2) * tau(2,c) & ! directions
                                    + Ji(i,j,k,e,3,3) * tau(3,c) ) ! transposed

                C0(i,j,c)   = tmp * ( sm(i,j,f,e,c)               &
                                    - sp(i,j,f,e,c)               &
                                    - 2 * cip * ( n1 * tau(1,c)   &
                                                + n2 * tau(2,c)   &
                                                + n3 * tau(3,c) ) )

              end do

            end do
            end do

            select case(P)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
#ifdef __LIBXSMM__
              call LIBXSMM_XMMCall(xmm, Ds_t, Ct, Ds_Ct)
#else
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
#endif
            end select

            do c = 1, 3
              do j = 0, P
              do i = 0, P

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
            do k = 0, P
            do j = 0, P
            do i = 0, P
              F_d(i,j,k,e,c) = F_d(i,j,k,e,c)           &
                             - ( Ds(0,i) * Cn(j,k,1,c)  &
                               + Ds(P,i) * Cn(j,k,2,c)  &
                               + Ds(0,j) * Cn(i,k,3,c)  &
                               + Ds(P,j) * Cn(i,k,4,c)  &
                               + Ds(0,k) * Cn(i,j,5,c)  &
                               + Ds(P,k) * Cn(i,j,6,c)  &
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
  !> Compose exterior traces at element boundary

  subroutine GetExteriorTraces(P, element, struct, bc, e, f, n, vm, sm, vp, sp)

    integer,               intent(in)    :: P       !< polynomial order
    class(MeshElement_3D), intent(in)    :: element !< element
    logical,               intent(in)    :: struct  !< F/T if un/structured
    character,             intent(in)    :: bc(:)   !< boundary conditions
    integer,               intent(in)    :: e       !< element ID
    integer,               intent(in)    :: f       !< element face
    real(RNP), contiguous, intent(in)    :: n (0:,0:,:,:,:) !< normals
    real(RNP), contiguous, intent(in)    :: vm(0:,0:,:,:,:) !< v⁻
    real(RNP), contiguous, intent(in)    :: sm(0:,0:,:,:,:) !< s⁻
    real(RNP), contiguous, intent(inout) :: vp(0:,0:,:,:,:) !< vb → v⁺
    real(RNP), contiguous, intent(inout) :: sp(0:,0:,:,:,:) !< sb → s⁺

    ! local variables ..........................................................

    real(RNP) :: n1, n2, n3, vm_n, vp_n
    integer   :: c, i, j, k, l, m

    i  = element % face(f) % i_neighbor

    if (i > 0) then

      ! copy fluxes from neighbor ..............................................

      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      if (struct) then
        do c = 1, 3
          do k = 0, P
          do j = 0, P
            vp(j,k,f,e,c) = vm(j,k,m,l,c)
            sp(j,k,f,e,c) = sm(j,k,m,l,c)
          end do
          end do
        end do
      else
        do c = 1, 3
          call element % AlignFromNeighborFace(f, i, vm(:,:,m,l,c), vp(:,:,f,e,c))
          call element % AlignFromNeighborFace(f, i, sm(:,:,m,l,c), sp(:,:,f,e,c))
        end do
      end if

    else

      ! apply boundary conditions ..............................................

      select case(bc(element % face(f) % boundary))

      case('D') ! Dirichlet
        do c = 1, 3
          do k = 0, P
          do j = 0, P
            vp(j,k,f,e,c) = -vm(j,k,f,e,c) + 2 * vp(j,k,f,e,c)
            sp(j,k,f,e,c) = -sm(j,k,f,e,c)
          end do
          end do
        end do

      case('S') ! free slip
        do k =0, P
        do j =0, P

          n1 = n(j,k,f,e,1)
          n2 = n(j,k,f,e,2)
          n3 = n(j,k,f,e,3)

          vm_n = n1 * vm(j,k,f,e,1) + n2 * vm(j,k,f,e,2) + n3 * vm(j,k,f,e,3)
          vp_n = n1 * vp(j,k,f,e,1) + n2 * vp(j,k,f,e,2) + n3 * vp(j,k,f,e,3)

          vp(j,k,f,e,1) =  vm(j,k,f,e,1) + 2 * n1 * (vp_n - vm_n)
          vp(j,k,f,e,2) =  vm(j,k,f,e,2) + 2 * n2 * (vp_n - vm_n)
          vp(j,k,f,e,3) =  vm(j,k,f,e,3) + 2 * n3 * (vp_n - vm_n)

          sp(j,k,f,e,1) = -sm(j,k,f,e,1)
          sp(j,k,f,e,2) = -sm(j,k,f,e,2)
          sp(j,k,f,e,3) = -sm(j,k,f,e,3)

        end do
        end do

      case default
        do c = 1, 3
          do k =0, P
          do j =0, P
            vp(j,k,f,e,c) =  vm(j,k,f,e,c)
            sp(j,k,f,e,c) = -sm(j,k,f,e,c)
          end do
          end do
        end do

      end select

    end if

  end subroutine GetExteriorTraces

  !-----------------------------------------------------------------------------
  !> Computation of tangential derivatives

  pure subroutine TangentialDerivatives(P, D, C, DC)
    integer,   intent(in)  :: P                !< polynomial order
    real(RNP), intent(in)  :: D  (0:P,0:P)     !< diff matrix
    real(RNP), intent(in)  :: C  (0:P,0:P,2,3) !< tangential contributions
    real(RNP), intent(out) :: DC (0:P,0:P,2,3) !< derivatives of C

    integer :: i, j, l, m

    do l = 1, 3
      do j = 0, P
      do i = 0, P
        DC(i,j,1,l) = 0
        DC(i,j,2,l) = 0
        do m = 0, P
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

  subroutine DispatchKernel_RSP(kernel, P, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: P
    integer(LIBXSMM_BLASINT_KIND), target, optional, intent(in) :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_FLOAT), target :: alpha = 1, beta = 0

    m = P + 1
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

  subroutine DispatchKernel_RDP(kernel, P, lda, ldb, ldc, flags, prefetch)

    type(C_FunPtr), intent(out) :: kernel
    integer, intent(in) :: P
    integer(LIBXSMM_BLASINT_KIND), intent(in), optional, target :: lda, ldb, ldc
    integer(C_INT), target, optional, intent(in) :: flags
    integer(C_INT), target, optional, intent(in) :: prefetch

    integer(LIBXSMM_BLASINT_KIND) :: m, n, k
    real(C_DOUBLE), target :: alpha = 1, beta = 0

    m = P + 1
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

end submodule MP_GetDiffusionTerm_DC
