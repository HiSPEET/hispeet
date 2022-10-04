!> summary:  3D DG elliptic operator: application with deformed mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/10/6
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Apply_DC
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

    module subroutine AddFluxes_DC_Gen(mesh, eop, a, Ji_n, nu, tr, r, f)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP),                       intent(in)    :: nu
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DC_Gen

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: LIBXSMM version

    module subroutine AddFluxes_DC_XSMM(mesh, eop, a, Ji_n, nu, tr, r, f)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP),                       intent(in)    :: nu
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DC_XSMM

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Application with irregular (deformed) mesh and constant diffusivity

  module subroutine Apply_DC(this, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    real(RNP),                       intent(in)  :: lambda     !< λ
    real(RNP),                       intent(in)  :: nu         !< ν
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: tr_buf

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of u and q_n

    integer :: po, ne, ng, np

    associate( mesh => this % sem % mesh           &
             , Jd   => this % sem % metrics % Jd   &
             , G    => this % sem % metrics % G    &
             , a    => this % sem % metrics % a    &
             , Ji_n => this % sem % metrics % Ji_n &
             , eop  => this % eop                  )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,2))
      tr_buf = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Elliptic( eop%w, eop%D, Jd, G & !  As r may contain ghost entries
                       , lambda, nu, u       & !  it is explicitly restricted to
                       , r(:,:,:,:ne), Ji_n  & !  local elements since no proper
                       , ub = tr(:,:,:,:,1)  & !  bound checking is performed in
                       , qb = tr(:,:,:,:,2)  ) !  TPO_Elliptic.

      ! transfer traces and apply boundary conditions ..........................

      call tr_buf % Transfer(mesh, tr, tag=1000)

      call EnforceBoundaryConditions( this, bv              &
                                    , jmp_u = tr(:,:,:,:,1) &
                                    , avg_q = tr(:,:,:,:,2) )

      call tr_buf % Merge(tr)

      ! add fluxes .............................................................

      select case(po)
      case(1:2)
        call AddFluxes_DC_Gen  (mesh, eop, a, Ji_n, nu, tr, r, f)
      case default
#ifdef __LIBXSMM__
        call AddFluxes_DC_XSMM (mesh, eop, a, Ji_n, nu, tr, r, f)
#else
        call AddFluxes_DC_Gen  (mesh, eop, a, Ji_n, nu, tr, r, f)
#endif
      end select

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_DC

  !=============================================================================
  ! Shared procedures

  !-----------------------------------------------------------------------------
  !> Compose element-boundary fluxes from flux traces

  subroutine GetBoundaryFluxes(element, struct, e, f, tr, jmp_u, avg_q)

    class(MeshElement_3D), intent(in) :: element
    logical,   intent(in)  :: struct        !< F/T for un/structured mesh
    integer,   intent(in)  :: e             !< element ID
    integer,   intent(in)  :: f             !< element face
    real(RNP), intent(in)  :: tr(:,:,:,:,:) !< traces of u, q_n
    real(RNP), intent(out) :: jmp_u(:,:)    !< normal jump n⋅[u]
    real(RNP), intent(out) :: avg_q(:,:)    !< average normal flux n⋅{q}

    contiguous :: tr, jmp_u, avg_q

    integer :: i, l, m

    i = element % face(f) % i_neighbor
    if (i > 0) then
      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      if (struct) then
        jmp_u = (tr(:,:,f,e,1) - tr(:,:,m,l,1))
        avg_q = (tr(:,:,f,e,2) - tr(:,:,m,l,2)) * HALF
      else
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,1), jmp_u)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,2), avg_q)
        jmp_u = (tr(:,:,f,e,1) - jmp_u)
        avg_q = (tr(:,:,f,e,2) - avg_q) * HALF
      end if
    else
      jmp_u = tr(:,:,f,e,1)
      avg_q = tr(:,:,f,e,2)
   end if

  end subroutine GetBoundaryFluxes

  !=============================================================================

end submodule MP_Apply_DC
