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

  module subroutine Stokes_V_Cycle(this, tau, mu, nu, bv, f, u, n_cyc, l_top)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), optional, intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution
    integer,  optional, intent(in) :: n_cyc
      !< number of cycles    [1]
    integer,  optional, intent(in) :: l_top
      !< top level        [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: r, v
    real(RNP), contiguous, pointer, save :: mu_l(:,:,:,:)
    real(RNP), contiguous, pointer, save :: nu_l(:,:,:,:)
    integer, save :: l_top_, n_cyc_

    integer :: l, m

    associate( problem => this % problem       &
             , sem     => this % ml_op_u % sem )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      call r % Init(this%ml_op_u, nc = problem%nc)
      call v % Init(this%ml_op_u, nc = problem%nc)

      mu_l => null()
      nu_l => null()

      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, n_cyc_

        V_DOWN: do l = l_top_, 2, -1

          associate( ins_l  => this % ins_op(l)      &
                   , mesh_l => sem(l  ) % mesh       &
                   , mesh_p => sem(l-1) % mesh       &
                   , bv_l   => bv % level(l  ) % var &
                   , f_l    => f  % level(l  ) % val &
                   , u_l    => u  % level(l  ) % val &
                   , r_l    => r  % level(l  ) % val )

            ! pre-smoothing and residual computation ...........................

            !$omp master
            if (present(mu)) mu_l => mu % level(l) % val(:,:,:,:,1)
            if (present(nu)) nu_l => nu % level(l) % val(:,:,:,:,1)
            !$omp end master
            !!omp barrier after initialization in StokesProjection/FGMRES

            call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)
            call ins_l % GetStokesResidual(tau, f_l, bv_l, mu_l, nu_l, u_l, r_l)

            ! restriction ......................................................
            ! parent RHS .......................................................

          end associate

        end do V_DOWN

      ! coarse grid solver ...............................................
        ! prolongation .................................................
        ! post-smoothing .................................................
      !?termination check ....................................................

      end do V_OUTER

      !$omp master
      deallocate(r, v)
      !$omp end master

    end associate

  end subroutine Stokes_V_Cycle

  !=============================================================================

end submodule
