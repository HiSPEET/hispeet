!> summary:  Application of essential boundary conditions
!> author:   Joerg Stiller
!> date:     2022/12/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyEssentialBC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of essential boundary conditions to velocity trace variables
  !>
  !> Construction of the outer boundary traces of velocity using inner traces
  !> and boundary values. In the absence of the latter, homogeneous conditions
  !> are assumed.

  module subroutine ApplyEssentialBC(this, bv_u, vm, vp)
    class(INS_Operator_3D), intent(in) :: this
    !< Navier-Stokes operator
    class(BoundaryVariable_3D), optional, intent(in) :: bv_u(:)
    !< boundary values depending on BC type specified in `this % bc_v`
    real(RNP), contiguous, intent(in) :: vm(:,:,:,:,:)
    !< inner velocity traces, `vm(np,np,6,ne,3) = v⁻`
    real(RNP), contiguous, intent(inout) :: vp(:,:,:,:,:)
    !< outer velocity traces, `vp(np,np,6,ne,3) = v⁺` on boundary faces,
    !! values on interior faces remain unchanged

    integer :: b, e, f, m

    do b = 1, this % mesh % n_bound
      associate(boundary => this % mesh % boundary(b))

        select case(this % bc_v(b))

        case('D') ! Dirichlet: v⁺ = 2vᵇ - v⁻

          !$omp do
          do f = 1, boundary % n_face
            e = boundary % face(f) % element_id
            m = boundary % face(f) % element_face
            if (present(bv_u)) then
              associate(vb => bv_u(b) % val)
                vp(:,:,m,e,1:3) = 2 * vb(:,:,f,1:3) - vm(:,:,m,e,1:3)
              end associate
            else
              vp(:,:,m,e,1:3) = -vm(:,:,m,e,1:3)
            end if
          end do

        case('O') ! Outflow: v⁺ = v⁻

          !$omp do
          do f = 1, boundary % n_face
            e = boundary % face(f) % element_id
            m = boundary % face(f) % element_face
            vp(:,:,m,e,1:3) = vm(:,:,m,e,1:3)
          end do

        end select

      end associate
    end do

  end subroutine ApplyEssentialBC

  !=============================================================================

end submodule MP_ApplyEssentialBC
