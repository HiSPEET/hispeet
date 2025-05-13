!> summary:  Stokes V-cycle for Stokes part in semi-implicit INS solvers
!> author:   Joerg Stiller
!> date:     2025/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_Stokes_V_Cycle
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Performs one or more FAS-MG V-cycles for the Stokes part

  module subroutine Stokes_V_Cycle(this, tau, mu, nu, f, u, n_cyc, l_top, r_2)
    class(ML_INS_Operator_3D), intent(in)    :: this
    real(RNP),                 intent(in)    :: tau   !< effective time step
    class(ML_MeshVariable_3D), intent(in)    :: mu    !< bulk viscosity
    class(ML_MeshVariable_3D), intent(in)    :: nu    !< shear viscosity
    class(ML_MeshVariable_3D), intent(in)    :: f     !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u     !< solution
    integer,         optional, intent(in)    :: n_cyc !< number of cycles    [1]
    integer,         optional, intent(in)    :: l_top !< top level        [auto]
    real(RNP),       optional, intent(out)   :: r_2   !< Euclidean residual norm

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: r, v
    integer, save :: l_top_, n_cyc_

    integer :: l, m

    associate( problem => this % ml_ins % problem            &
             , sem     => this % ml_ins % ml_op_u % sem      &
             )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      if (present(n_cyc)) then
        n_cyc_ = n_cyc
      else
        n_cyc_ = 1
      end if
      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if
      allocate(r, v)
      call r % Init(this%ml_op, nc = problem%nc)
      call v % Init(this%ml_op, nc = problem%nc)
      !$omp end master

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, n_cyc_

        V_DOWN: do l = l_top_, 2, -1

          associate( mesh_l => sem(l  ) % mesh                  &
                   , mesh_p => sem(l-1) % mesh                  &
                   )

            ! pre-smoothing and residual computation ...........................

            call ins % StokesSolver(tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)

        ! restriction ....................................................
        ! parent RHS .....................................................
      ! coarse grid solver ...............................................
        ! prolongation .................................................
        ! post-smoothing .................................................
      !?termination check ....................................................

    end associate

  end subroutine Stokes_V_Cycle

  !=============================================================================

end submodule
