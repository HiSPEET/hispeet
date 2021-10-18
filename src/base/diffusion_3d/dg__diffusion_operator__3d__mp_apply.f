!> summary:  3D DG diffusion operator: application
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D) MP_Apply
  use Mesh__3D
  use Mesh_Element__3D
  implicit none

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Application with regular mesh and constant diffusivity

    module subroutine Apply_RC(this, u, v, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      real(RNP), intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), intent(out) :: v(:,:,:,:) !< result
      real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS
    end subroutine Apply_RC

    !---------------------------------------------------------------------------
    !> Application with regular mesh and variable diffusivity

    module subroutine Apply_RV(this, u, v, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      real(RNP), intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), intent(out) :: v(:,:,:,:) !< result
      real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS
    end subroutine Apply_RV

    !---------------------------------------------------------------------------
    !> Application with irregular (deformed) mesh and constant diffusivity

    module subroutine Apply_DC(this, u, v, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      real(RNP), intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), intent(out) :: v(:,:,:,:) !< result
      real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS
    end subroutine Apply_DC

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator

  module subroutine Apply(this, u, v, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result
    real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS

    if (this % sem % mesh % regular) then
      if (allocated(this % nu_pv)) then
        ! regular variable
        call Apply_RV(this, u, v, f)
      else
        ! regular constant
        call Apply_RC(this, u, v, f)
      end if
    else
      ! deformed constant
      call Apply_DC(this, u, v, f)
    end if

  end subroutine Apply

  !=============================================================================
  ! Shared procedures

  !-----------------------------------------------------------------------------
  !> Compose element-boundary fluxes from flux traces for homogeneous BC

  subroutine GetElementBoundaryFluxes(element, e, f, tr_u, tr_qn, Ju, Aq)
    class(MeshElement_3D), intent(in)  :: element        !< element
    integer,               intent(in)  :: e              !< element ID
    integer,               intent(in)  :: f              !< element face
    real(RNP), contiguous, intent(in)  :: tr_u (:,:,:,:) !< u   @ element faces
    real(RNP), contiguous, intent(in)  :: tr_qn(:,:,:,:) !< q_n @ element faces
    real(RNP), contiguous, intent(out) :: Ju(:,:) !< normal jump n⋅[u]  w/o BC
    real(RNP), contiguous, intent(out) :: Aq(:,:) !< average flux n⋅{q} w/o BC

    integer :: i, l, m

    i = element % face(f) % i_neighbor
    if (i > 0) then
      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      Ju = (tr_u (:,:,f,e) - tr_u (:,:,m,l))
      Aq = (tr_qn(:,:,f,e) - tr_qn(:,:,m,l)) * HALF
    else
      Ju = tr_u (:,:,f,e)
      Aq = tr_qn(:,:,f,e)
   end if

  end subroutine GetElementBoundaryFluxes

  !-----------------------------------------------------------------------------
  !> Modify boundary traces to yield correct contribution to the operator
  !>
  !> On Dirichlet boundaries (bc = 'D'):
  !>
  !>   – interior solution contributes twice: [u]ᵢ = 2u⁻
  !>   - du/dn is extrapolated from interior: {q}ᵢ = qᵢ⁻
  !>
  !> and on Neumann boundaries (bc = 'N'):
  !>
  !>   – solution is extrapolated from interior: [u]ᵢ = 0
  !>   - du/dn is reflected from interior:       {q}ᵢ = 0
  !>
  !> where `i` is the coordinate direction parallel to the face normal.

  subroutine ApplyBoundaryConditions(mesh, bc, tr_u, tr_qn)
    class(Mesh_3D), intent(in)    :: mesh           !< mesh partition
    character,      intent(in)    :: bc(:)          !< boundary conditions
    real(RNP),      intent(inout) :: tr_u (:,:,:,:) !< trace of u
    real(RNP),      intent(inout) :: tr_qn(:,:,:,:) !< trace of ν du/dn

    integer :: b, e, f, l

    do b = 1, size(bc)
      associate(boundary_face => mesh % boundary(b) % face)

        select case(bc(b))

        case('D') ! Dirichlet
          !$omp do
          do l = 1, size(boundary_face)
            e = boundary_face(l) % mesh_element % id    ! mesh element
            f = boundary_face(l) % mesh_element % face  ! element face
            tr_u (:,:,f,e) = 2 * tr_u (:,:,f,e)
          end do
          !$omp end do nowait

        case('N') ! Neumann
          ! Neumann: du/dn does not contribute, [u] = 0 due to extrapolation
          !$omp do
          do l = 1, size(boundary_face)
            e = boundary_face(l) % mesh_element % id    ! mesh element
            f = boundary_face(l) % mesh_element % face  ! element face
            tr_u (:,:,f,e) = 0
            tr_qn(:,:,f,e) = 0
          end do
          !$omp end do nowait

        end select
      end associate
    end do

    !$omp barrier

  end subroutine ApplyBoundaryConditions

  !=============================================================================

end submodule MP_Apply
