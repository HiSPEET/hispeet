!> summary:  3D DG elliptic operator: evaluation with deformed mesh and
!>           variable isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2024/08/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Eval_DV
  use TPO__Elliptic__3D
  use Mesh_Metrics__3D
  use Mesh_Element__3D
  use Mesh__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: generic version

    module subroutine AddFluxes_DV_Gen(mesh, eop, a, Ji_n, tr, r, f)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DV_Gen

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: LIBXSMM version

    module subroutine AddFluxes_DV_XSMM(mesh, eop, a, Ji_n, tr, r, f)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DV_XSMM

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Evaluation with deformed mesh and variable isotropic diffusivity
  !>
  !> Usage
  !>   1) `f` and `bv` given:  computation of the residual, `r = f - Au`
  !>   2) `f` and `bv` absent: evaluation of the homogeneous operator, `r = Au`

  module subroutine Eval_DV(this, bc, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    character,                       intent(in)  :: bc(:)       !< {P,D,N}
    real(RNP),                       intent(in)  :: lambda      !< λ
    real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:)  !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:)  !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:)  !< RHS
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of ν, u and q_n

    integer :: po, ne, ng, np

    if (this % sem % mesh % part < 0) return

    associate( eop  => this % eop                  &
             , sem  => this % sem                  &
             , mesh => this % sem % mesh           &
             , a    => this % sem % metrics % a    &
             , Ji_n => this % sem % metrics % Ji_n )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,3))
      buf_tr = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Elliptic( eop, sem            & ! As r may contain ghost entries
                       , lambda, nu, u       & ! it is explicitly restricted to
                       , v   = r(:,:,:,:ne)  & ! local elements since no proper
                       , nub = tr(:,:,:,:,1) & ! bound checking is performed in
                       , ub  = tr(:,:,:,:,2) & ! bound checking is performed in
                       , qb  = tr(:,:,:,:,3) ) ! TPO_Elliptic.

      ! transfer traces and apply boundary conditions ..........................

      call buf_tr % Transfer(mesh, tr, tag=1000)

      call EnforceBoundaryConditions( this, bc, bv          &
                                    , jmp_u = tr(:,:,:,:,2) &
                                    , avg_q = tr(:,:,:,:,3) )

      call buf_tr % Merge(tr)

      ! add fluxes .............................................................

      select case(po)
      case(1:2)
        call AddFluxes_DV_Gen  (mesh, eop, a, Ji_n, tr, r, f)
      case default
#ifdef __LIBXSMM__
        call AddFluxes_DV_XSMM (mesh, eop, a, Ji_n, tr, r, f)
#else
        call AddFluxes_DV_Gen  (mesh, eop, a, Ji_n, tr, r, f)
#endif
      end select

      ! frozen and ghost elements are set to zero ..............................

      call SetArray(r(:,:,:,mesh%n_elem_active+1:), ZERO)

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, buf_tr)
      !$omp end master

    end associate

  end subroutine Eval_DV

  !=============================================================================
  ! Shared procedures

  !-----------------------------------------------------------------------------
  !> Compose element-boundary fluxes from flux traces

  subroutine GetBoundaryFluxes(element, struct, e, f, tr, nu_max, jmp_u, avg_q)

    class(MeshElement_3D), intent(in) :: element
    logical,   intent(in)  :: struct        !< F/T for un/structured mesh
    integer,   intent(in)  :: e             !< element ID
    integer,   intent(in)  :: f             !< element face
    real(RNP), intent(in)  :: tr(:,:,:,:,:) !< traces of u, q_n
    real(RNP), intent(out) :: nu_max(:,:)   !< max(ν⁻,ν⁺)
    real(RNP), intent(out) :: jmp_u(:,:)    !< normal jump n⋅[u]
    real(RNP), intent(out) :: avg_q(:,:)    !< average normal flux n⋅{q}

    contiguous :: tr, jmp_u, avg_q

    integer :: i, l, m

    i = element % face(f) % i_neighbor
    if (i > 0) then
      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      if (struct) then
        nu_max = max(tr(:,:,f,e,1), tr(:,:,m,l,1))
        jmp_u = (tr(:,:,f,e,2) - tr(:,:,m,l,2))
        avg_q = (tr(:,:,f,e,3) - tr(:,:,m,l,3)) * HALF
      else
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,1), nu_max)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,2), jmp_u)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,3), avg_q)
        jmp_u = (tr(:,:,f,e,2) - jmp_u)
        avg_q = (tr(:,:,f,e,3) - avg_q) * HALF
      end if
    else
      nu_max = tr(:,:,f,e,1)
      jmp_u  = tr(:,:,f,e,2)
      avg_q  = tr(:,:,f,e,3)
   end if

  end subroutine GetBoundaryFluxes

  !=============================================================================

end submodule MP_Eval_DV
