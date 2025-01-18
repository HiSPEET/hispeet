!> summary:  Outflow pressure penalty opposing reverse flow
!> author:   Joerg Stiller
!> date:     2023/10/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetBackflowPenalty

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Backflow pressure penalty

  module subroutine GetBackflowPenalty(this, problem, b, v, dp)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    class(INS_Problem_3D), intent(in) :: problem
    !< incompressible Navier-Stokes problem

    integer, intent(in) :: b
    !< boundary ID

    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)
    !< velocity

    real(RNP), contiguous, intent(out) :: dp(0:,0:,:)
    !< pressure correction opposing reverse flow at outlets, dp = -E(v,n)

    ! local variables ..........................................................

    real(RNP), parameter :: eps = epsilon(ONE)
    real(RNP) :: c_vn, v1, v2, v3, vn, vv
    integer   :: e, f, i, j, m

    ! scaling factor
    c_vn = 1 / max(eps, this % delta_outflow * problem % v_ref)

    associate( boundary => this % sem_u % mesh % boundary(b) &
             , n        => this % sem_u % metrics % n        &
             , po       => this % eop_u % po                 )

      !$omp do
      do f = 1, boundary % n_face
        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

        do j = 0, po
        do i = 0, po

          ! velocity at current position
          select case(m)
          case(1:2)
            v1 = v((m-1)*po,i,j,e,1)
            v2 = v((m-1)*po,i,j,e,2)
            v3 = v((m-1)*po,i,j,e,3)
          case(3:4)
            v1 = v(i,(m-1)*po,j,e,1)
            v2 = v(i,(m-1)*po,j,e,2)
            v3 = v(i,(m-1)*po,j,e,3)
          case(5:6)
            v1 = v(i,j,(m-1)*po,e,1)
            v2 = v(i,j,(m-1)*po,e,2)
            v3 = v(i,j,(m-1)*po,e,3)
          end select

          ! normal and squared velocity
          vn = v1 * n(i,j,m,e,1) + v2 * n(i,j,m,e,2) + v3 * n(i,j,m,e,3)
          vv = v1 * v1 + v2 * v2 + v3 * v3

          ! pressure correction, ∆p = -E(v,n) ≤ 0
          if (vn < 0) then
            dp(i,j,f) = -HALF * vv * tanh(c_vn * abs(vn))
          else
            dp(i,j,f) = 0
          end if

        end do
        end do

      end do

    end associate

  end subroutine GetBackflowPenalty

  !=============================================================================

end submodule MP_GetBackflowPenalty
