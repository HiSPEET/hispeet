!> summary:  3D DG elliptic operator: application
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Apply
  use Mesh__3D
  use Mesh_Element__3D
  implicit none

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Application with regular mesh and constant diffusivity

    module subroutine Apply_RC(this, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
      !< boundary values
    end subroutine Apply_RC

    !---------------------------------------------------------------------------
    !> Application with regular mesh and variable diffusivity

    module subroutine Apply_RV(this, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
      !< boundary values
    end subroutine Apply_RV

    !---------------------------------------------------------------------------
    !> Application with irregular (deformed) mesh and constant diffusivity

    module subroutine Apply_DC(this, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
      !< boundary values
    end subroutine Apply_DC

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator

  module subroutine Apply(this, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    if (this % sem % mesh % regular) then
      if (allocated(this % nu_pv)) then
        ! regular variable
        call Apply_RV(this, u, r, f, bv)
      else
        ! regular constant
        call Apply_RC(this, u, r, f, bv)
      end if
    else
      ! deformed constant
      call Apply_DC(this, u, r, f, bv)
    end if

  end subroutine Apply

  !=============================================================================
  ! Shared procedures

  !-----------------------------------------------------------------------------
  !> Weak enforcement of boundary conditions
  !>
  !> Input:
  !>
  !>   - `bv`     Dirichlet or Neumann boundary values in component 1
  !>   - `jmp_u`  u⁻    on boundary element faces
  !>   - `avg_q`  n⋅q⁻  on boundary element faces
  !>
  !> Interior and ghost face entries `jmp_u` and `avg_q` will be ignored and
  !> remain unchanged.
  !>
  !> On output, the boundary conditions are applied as follows:
  !>
  !>   - Dirichlet boundaries (bc = 'D'):
  !>
  !>          jmp_u  ←  n⋅[u]  =  2(u⁻ - u_b)
  !>          avg_q  ←  n⋅{q}  =  n⋅q⁻            (unchanged)
  !>
  !>   - Neumann boundaries (bc = 'N'):
  !>
  !>          jmp_u  ←  n⋅[u]  =  0
  !>          avg_q  ←  n⋅{q}  =  q_b

  subroutine EnforceBoundaryConditions(diffusion_op, bv, jmp_u, avg_q)
    class(DG_EllipticOperator_3D), intent(in) :: diffusion_op
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(inout) :: jmp_u(:,:,:,:) !< trace of u
    real(RNP), contiguous, intent(inout) :: avg_q(:,:,:,:) !< trace of ν du/dn

    logical :: has_bv
    integer :: b, e, f, l

    has_bv = present(bv)

    associate(boundary => diffusion_op % sem % mesh % boundary)

      do b = 1, size(boundary)


        select case(diffusion_op % bc(b))

        case('D')  ! avg_q remains unchanged !

          if (has_bv) then
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = 2 * (jmp_u (:,:,f,e) - bv(b) % val(:,:,l,1))
            end do
            !$omp end do nowait
          else
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = 2 * jmp_u (:,:,f,e)
            end do
            !$omp end do nowait
          end if

        case('N')

          if (has_bv) then
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = ZERO
              avg_q(:,:,f,e) = bv(b) % val(:,:,l,1)
            end do
            !$omp end do nowait
          else
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = ZERO
              avg_q(:,:,f,e) = ZERO
            end do
            !$omp end do nowait
          end if

        end select

      end do
     !$omp barrier

    end associate

  end subroutine EnforceBoundaryConditions

  !=============================================================================

end submodule MP_Apply
