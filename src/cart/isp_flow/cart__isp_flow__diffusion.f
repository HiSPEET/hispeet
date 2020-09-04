!> summary:  ISP flow: implicit diffusion step
!> author:   Joerg Stiller
!> date:     2018/04/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: implicit diffusion step
!===============================================================================

module CART__ISP_Flow__Diffusion

  use Kind_Parameters, only: RNP
  use TPO__Diagonal_3d
  use ISP_Flow_Problem
  use CART__Boundary_Variable
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: DiffusionStep

contains

!-------------------------------------------------------------------------------
!>  Diffusion step

subroutine DiffusionStep(problem, flow_op, dt, f, u, w, nu, i_max)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem        !< flow problem
  class(FlowOperators),    intent(inout) :: flow_op        !< flow operators
  real(RNP),               intent(in)    :: dt             !< time step size
  real(RNP),               intent(in)    :: f(:,:,:,:,:)   !< sources
  real(RNP),               intent(inout) :: u(:,:,:,:,:)   !< solution
  real(RNP),               intent(inout) :: w(:,:,:,:,:)   !< workspace
  real(RNP),     optional, intent(in)    :: nu(:,:,:,:,:)  !< diffusivities
  integer,       optional, intent(in)    :: i_max          !< max num iterations

  ! local variables ............................................................

  real(RNP) :: g, r_2
  integer   :: ne, ni, np
  integer   :: c

  associate( mesh => flow_op % mesh   &
           , eop  => flow_op % eop_u  &
           , pmg  => flow_op % pmg_u  &
           , bc   => problem % bc     &
           , fc   => w(:,:,:,:,1)     )

    ! intialization ............................................................

    ! dimensions
    np = size(u, 1)
    ne = mesh % ne

    ! solve components .........................................................

    Components: do c = 1, size(u,5)

      if (c == 4) cycle ! skip pressure

      ! set problem
      if (present(nu)) then
        call pmg % SetProblem(1/dt, nu(:,:,:,:,c), bc(:,c))
      else
        if (allocated(problem%nu_svv_ref)) then
          call pmg % SetProblem(1/dt,                  &
                                problem%nu_ref(c),     &
                                problem%nu_svv_ref(c), &
                                bc(:,c)                )
        else
          call pmg % SetProblem(1/dt, problem%nu_ref(c), bc(:,c))
        end if
      end if

      ! project sources
      g = product(mesh%dx) / (8 * dt)
      call TPO_Diagonal(g, eop%w, f(:,:,:,:,c), fc)

      ! add boundary contributions
      call pmg % BcToRHS(flow_op%bv_u, c, fc)

      ! solve
      call pmg % MG_CG_Solver(u(:,:,:,:,c), fc, ni, r_2, i_max)

      ! monitoring
      !$omp single
      if (flow_op % control % monitor > 0 .and. mesh % part == 0) then
        print '(4X,A,I0,A,I4,A,ES9.2)', &
              'diffusion u[',c,']: ni =', ni, ', ‖r‖ =',r_2
      end if
      !$omp end single

    end do Components

  end associate

end subroutine DiffusionStep

!===============================================================================

end module CART__ISP_Flow__Diffusion
