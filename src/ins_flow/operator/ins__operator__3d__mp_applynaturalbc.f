!> summary:  Application of natural boundary conditions
!> author:   Joerg Stiller
!> date:     2023/10/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyNaturalBC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of natural boundary conditions to velocity trace variables
  !>
  !> The outer viscous flux vector traces s⁺ = n⁺⋅τ⁺ are constructed using the
  !> inner fluxes s⁻ and the boundary values defined in `bv_s`. If the latter
  !> are not present, homogeneous conditions are assumed.

  module subroutine ApplyNaturalBC(this, bv_s, sm, sp)
    class(INS_Operator_3D), intent(in) :: this
    !< Navier-Stokes operator
    class(BoundaryVariable_3D), intent(in) :: bv_s(:)
    !< viscous boundary fluxes depending on BC type
    !!   - Γᴰ :  not used
    !!   - Γᴼ :  [ s₁, s₂, s₃ ]   viscous fluxes
    real(RNP), contiguous, intent(in) :: sm(:,:,:,:,:)
    !< inner viscous flux vector, `sm(np,np,6,ne,3) = s⁻ = n⁻⋅τ⁻`
    real(RNP), contiguous, intent(inout) :: sp(:,:,:,:,:)
    !< outer viscous flux vector, `sp(np,np,6,ne,3) = s⁺ = n⁺⋅τ⁺` on boundary
    !! faces, values on interior faces remain unchanged

    integer :: b, e, f, m

    do b = 1, this % mesh % n_bound
      associate(boundary => this % mesh % boundary(b))

        select case(this % problem % bc_v(b))

        case('D') ! Dirichlet: s⁺ = -s⁻

          !$omp do
          do f = 1, boundary % n_face
            e = boundary % face(f) % element_id
            m = boundary % face(f) % element_face
            sp(:,:,m,e,1:3) = -sm(:,:,m,e,1:3)
          end do

        case('O') ! Outflow: s⁺ = s⁻ - 2sᵇ

          associate(sb => bv_s(b) % val)
            !$omp do
            do f = 1, boundary % n_face
              e = boundary % face(f) % element_id
              m = boundary % face(f) % element_face
              sp(:,:,m,e,1:3) = sm(:,:,m,e,1:3) - 2 * sb(:,:,f,1:3)
            end do
          end associate

        end select

      end associate
    end do

  end subroutine ApplyNaturalBC

  !=============================================================================

end submodule MP_ApplyNaturalBC
