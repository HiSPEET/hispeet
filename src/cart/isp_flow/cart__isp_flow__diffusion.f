!> summary:  ISP flow: implicit diffusion step
!> author:   Joerg Stiller
!> date:     2018/04/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: implicit diffusion step
!===============================================================================

module CART__ISP_Flow__Diffusion

  use Kind_Parameters, only: RNP
  use TPO_sDDD
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

subroutine DiffusionStep(problem, flow_op, dt, f, u, w, nu)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem        !< flow problem
  class(FlowOperators),    intent(inout) :: flow_op        !< flow operators
  real(RNP),               intent(in)    :: dt             !< time step size
  real(RNP),               intent(in)    :: f(:,:,:,:,:)   !< sources
  real(RNP),               intent(inout) :: u(:,:,:,:,:)   !< solution
  real(RNP),               intent(inout) :: w(:,:,:,:,:)   !< workspace
  real(RNP),     optional, intent(in)    :: nu(:,:,:,:,:)  !< diffusivities

  ! local variables ............................................................

  type(BoundaryVariable), allocatable, save :: bv_uc(:)

  real(RNP) :: g, r_2
  integer   :: l_top, ne, ni, np
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
    l_top = ubound(pmg%level, 1)

    ! workspace
    !$omp single
    allocate(bv_uc(mesh % n_boundary))
    !$omp end single

    ! solve components .........................................................

    Components: do c = 1, size(u,5)

      if (c == 4) cycle ! skip pressure

      ! set problem
      if (present(nu)) then
        call pmg % SetProblem(1/dt, nu(:,:,:,:,c), bc(:,c))
      else
        call pmg % SetProblem(1/dt, problem%nu_ref(c), bc(:,c))
      end if

      ! project sources
      g = product(mesh%dx) / 8
      call TPO_sDDD_Eval(np, ne, g, eop%w, f(:,:,:,:,c), fc)

      ! add boundary contributions
      call flow_op % bv_u % GetHandle(c, bv_uc)
      call pmg % level(l_top) % elliptic_op % BcToRHS(bv_uc, fc)

      ! solve
      call pmg % MG_CG_Solver(u(:,:,:,:,c), fc, ni, r_2)

      ! monitoring
      !$omp single
      if (flow_op % monitor_level > 0 .and. mesh % part == 0) then
        print '(4X,A,I0,A,I4,A,ES9.2)','diffusion u[',c,']: ni',ni,', ‖r‖ =',r_2
      end if
      !$omp end single

    end do Components

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(bv_uc)
    !$omp end master

  end associate

end subroutine DiffusionStep

!===============================================================================

end module CART__ISP_Flow__Diffusion
