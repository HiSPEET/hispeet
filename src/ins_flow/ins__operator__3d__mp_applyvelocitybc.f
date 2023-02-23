!> summary:  Application of velocity boundary conditions
!> author:   Joerg Stiller
!> date:     2022/12/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyVelocityBC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of velocity boundary conditions to trace variables
  !>
  !> The boundary variable `bv_v` is expected to contain the velocity boundary
  !> values in the first three components. Velocity values will be injected
  !> into `tr_v` and stresses into `tr_s`, if present.
  !>
  !> The traces are stored as element face variables defined as
  !>
  !>     tr_v(np,np,6,nl,3)
  !>     tr_s(np,np,6,nl,3)
  !>
  !> where `np` is the number of velocity element points per direction and
  !> `nl` is the number of elements, possibly including the ghosts.
  !>
  !> @note
  !> So far, only Dirichlet conditions ('D') are supported.

  module subroutine ApplyVelocityBC(this, bv_v, tr_v, tr_s)
    class(INS_Operator_3D), intent(in) :: this
    !< Navier-Stokes operator
    class(BoundaryVariable_3D), optional, intent(in) :: bv_v(:)
    !< boundary values
    real(RNP), contiguous, optional, intent(inout) :: tr_v(:,:,:,:,:)
    !< velocity on element faces ∊ ∂Ω: v⁻ → v⁺
    real(RNP), contiguous, optional, intent(inout) :: tr_s(:,:,:,:,:)
    !< stress vector on element faces ∊ ∂Ω: s⁻ → s⁺

    integer :: b, c, e, f, m

    do b = 1, this % mesh % n_bound
      associate(boundary => this % mesh % boundary(b))

        if (this % bc_v(b) == 'D') then

          ! Dirichlet conditions: v⁺ = 2vᵇ - v⁻, s⁺ = s⁻

          if (present(tr_v)) then
            if (present(bv_v)) then
              !$omp do
              do f = 1, boundary % n_face
                e = boundary % face(f) % element_id
                m = boundary % face(f) % element_face
                do c = 1, 3
                  tr_v(:,:,m,e,c) = 2 * bv_v(b) % val(:,:,f,c) - tr_v(:,:,m,e,c)
                end do
              end do
            else
              !$omp do
              do f = 1, boundary % n_face
                e = boundary % face(f) % element_id
                m = boundary % face(f) % element_face
                do c = 1, 3
                  tr_v(:,:,m,e,c) = -tr_v(:,:,m,e,c)
                end do
              end do
            end if
          end if

        end if

      end associate
    end do

    if (present(tr_s)) then
      return ! nothing to do yet
    end if

  end subroutine ApplyVelocityBC

  !=============================================================================

end submodule MP_ApplyVelocityBC
