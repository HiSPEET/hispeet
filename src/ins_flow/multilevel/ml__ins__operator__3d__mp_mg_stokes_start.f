!> summary:  Cascade and FMG start procedures for Stokes multigrid solver
!> author:   Joerg Stiller
!> date:     2025/05/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_MG_Stokes_Start
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Cascade and FMG start for the Stokes multigrid solver

  module subroutine MG_Stokes_Start(this, tau, mu, nu, bv, f_d0, f, u, n_cyc)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(in) :: f_d0
      !< unweighted approximate diffusion term
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< unweighted RHS: f = v₀/τ + f_c + f_s + ...
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution
    integer,  optional, intent(in) :: n_cyc
      !< number of V-cycles before advancing to next level  [0]

    character(len=:), allocatable :: log_prefix
    integer :: l, l_top

    !$omp master
    associate(mesh => this % ins_op(1) % mesh)
      if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1) then
        log_prefix = LoggingPrefix('MG_Stokes_Start', mesh%part, mesh%n_parts)
      end if
    end associate
    !$omp end master

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'start'
    end if

    l_top = size(u%level)

    do l = 1, l_top

      call this % ins_op(l) % StokesSolver( tau                            &
                                          , f    % level(l)%val            &
                                          , bv   % level(l)%var            &
                                          , mu   % level(l)%val(:,:,:,:,1) &
                                          , nu   % level(l)%val(:,:,:,:,1) &
                                          , u    % level(l)%val            &
                                          , f_d0 % level(l)%val            )

      if (l == l_top) exit

      if (present(n_cyc)) then
        if (n_cyc > 0) then
          call this % MG_Stokes_Cycle(tau, mu, nu, bv, f, u, n_cyc, l_top = l)
        end if
      end if

      call ParentToChildInterpolation_3D                   &
               ( parent = this % ml_op_u % sem(l  ) % mesh &
               , child  = this % ml_op_u % sem(l+1) % mesh &
               , iop    = this % ml_op_u % iop_cf_x(l)     &
               , v_p    = u % level(l  ) % val             &
               , v_c    = u % level(l+1) % val             )

    end do

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'exit'
    end if

  end subroutine MG_Stokes_Start

  !=============================================================================

end submodule MP_MG_Stokes_Start
