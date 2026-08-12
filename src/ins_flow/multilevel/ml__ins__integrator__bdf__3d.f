!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Multilevel IMEX BDF method for incompressible flows
!> author:   Joerg Stiller
!> date:     2025/05/12
!===============================================================================

module ML__INS__Integrator__BDF__3D
  use Kind_Parameters, only: RNP
  use Constants
  use Logging_Levels
  use XMPI
  use Array_Assignments
  use Boundary_Variable__3D
  use INS__Problem__3D
  use INS__Integrator__BDF__PrepStep__3D
  use Parent_To_Child_Interpolation__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__INS__Operator__3D
  use ML__INS__Integrator__3D
  implicit none
  private

  public :: ML_INS_Integrator_BDF_3D
  public :: ML_INS_Integrator_BDF_Options_3D

  !-----------------------------------------------------------------------------
  !> Type providing IMEX BDF solvers for 3D incompressible flows

  type, extends(ML_INS_Integrator_3D) :: ML_INS_Integrator_BDF_3D
    logical :: coupled !< switch for using coupled Stokes solver
    integer :: n_fmg   !< num cycles used in FMG start, 0 for cascade
    integer :: n_cyc   !< num cycles used by coupled Stokes solver
  contains
    procedure, non_overridable :: Init_ML_INS_Integrator_BDF_3D
    procedure :: TimeStep
  end type ML_INS_Integrator_BDF_3D

  ! constructor
  interface ML_INS_Integrator_BDF_3D
    module procedure New_ML_INS_Integrator_BDF_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing multilevel BDF options (none, so far)

  type, extends(ML_INS_IntegratorOptions_3D) :: ML_INS_Integrator_BDF_Options_3D
    logical :: coupled = .false. !< switch for using coupled Stokes solver
    integer :: n_fmg   = 0       !< num cycles used in FMG start, 0 for cascade
    integer :: n_cyc   = 1       !< num cycles used by coupled Stokes solver
  contains
    procedure :: Bcast => Bcast_ML_INS_Integrator_BDF_Options
  end type ML_INS_Integrator_BDF_Options_3D

contains

  !=============================================================================
  ! TBP of ML_INS_Integrator_BDF_Options_3D

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type ML_INS_Integrator_BDF_3D

  function New_ML_INS_Integrator_BDF_3D(problem, ml_ins, opt) result(this)
    class(INS_Problem_3D),                   intent(in) :: problem
    class(ML_INS_Operator_3D),               intent(in) :: ml_ins
    class(ML_INS_Integrator_BDF_Options_3D), intent(in) :: opt
    type(ML_INS_Integrator_BDF_3D) :: this

    call Init_ML_INS_Integrator_BDF_3D(this, problem, ml_ins, opt)

  end function New_ML_INS_Integrator_BDF_3D

  !-----------------------------------------------------------------------------
  !> Initialization of a ML_INS_Integrator_BDF_3D object

  subroutine Init_ML_INS_Integrator_BDF_3D(this, problem, ml_ins, opt)
    class(ML_INS_Integrator_BDF_3D),      intent(inout) :: this
    class(INS_Problem_3D),           target, intent(in) :: problem
    class(ML_INS_Operator_3D),       target, intent(in) :: ml_ins
    class(ML_INS_Integrator_BDF_Options_3D), intent(in) :: opt

    ! intialize parent type
    call this % Init_ML_INS_Integrator_3D(problem, ml_ins, opt)

    this % coupled = opt % coupled
    this % n_fmg   = opt % n_fmg
    this % n_cyc   = opt % n_cyc

  end subroutine Init_ML_INS_Integrator_BDF_3D

  !-----------------------------------------------------------------------------
  !> Execution of a multilevel BDF time step

  subroutine TimeStep(this, t, dt, u, first, last)
    class(ML_INS_Integrator_BDF_3D), intent(inout) :: this
    real(RNP),                 intent(inout) :: t     !< time t₀ → t
    real(RNP),                 intent(in)    :: dt    !< step size ∆t = t-t₀
    class(ML_MeshVariable_3D), intent(inout) :: u     !< u(x,t₀) → u(x,t)
    logical,         optional, intent(in)    :: first !< T for first step [F]
    logical,         optional, intent(in)    :: last  !< T for last  step [F]

    ! internal variables .......................................................

    class(ML_MeshVariable_3D), allocatable, save :: f   ! unwtd RHS
    class(ML_MeshVariable_3D), allocatable, save :: f_d ! unwtd diffusion term
    class(ML_MeshVariable_3D), allocatable, save :: mu  ! bulk diffusivity μ
    class(ML_MeshVariable_3D), allocatable, save :: nu  ! shear diffusivity ν

    ! saved at time t₁ = t₀-∆t
    class(ML_MeshVariable_3D), allocatable, save :: u1   ! flow variables
    class(ML_MeshVariable_3D), allocatable, save :: f_c1 ! unwtd convection term
    class(ML_MeshVariable_3D), allocatable, save :: f_d1 ! unwtd diffusion term

    class(ML_BoundaryVariable_3D), allocatable, save :: bv ! boundary values

    character(len=:), allocatable :: log_prefix
    real(RNP) :: tau
    integer   :: c, e, l

    associate( problem => this % problem &
             , ml_ins  => this % ml_ins  )

      ! initialization .........................................................

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      associate(mesh => ml_ins % ins_op(1)%mesh)
        if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1 ) then
          log_prefix = LoggingPrefix( 'TimeStep', mesh%part, mesh%n_parts)
        end if
      end associate

      if (allocated(log_prefix)) then
        print '(2A)', log_prefix, 'start'
      end if

      if (first) then
        if (allocated( f    )) deallocate( f    )
        if (allocated( f_d  )) deallocate( f_d  )
        if (allocated( mu   )) deallocate( mu   )
        if (allocated( nu   )) deallocate( nu   )
        if (allocated( u1   )) deallocate( u1   )
        if (allocated( f_c1 )) deallocate( f_c1 )
        if (allocated( f_d1 )) deallocate( f_d1 )
        if (allocated( bv   )) deallocate( bv   )
      end if
      if (.not. allocated(f)) then
        allocate(f, f_d, mu, nu, u1, f_c1, f_d1, bv)
        call f    % Init( ml_ins % ml_op_u, nc = problem % nc )
        call f_d  % Init( ml_ins % ml_op_u, nc = problem % nc )
        call mu   % Init( ml_ins % ml_op_u, nc = 1            )
        call nu   % Init( ml_ins % ml_op_u, nc = 1            )
        call u1   % Init( ml_ins % ml_op_u, nc = problem % nc )
        call f_c1 % Init( ml_ins % ml_op_u, nc = problem % nc )
        call f_d1 % Init( ml_ins % ml_op_u, nc = problem % nc )
        call bv   % Init( ml_ins % ml_op_u, nc = problem % nc )
      end if
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!omp barrier    !! not required because of barrier in BDF_PrepStep_3D

      ! preparation and extrapolation ..........................................

      do l = 1, size(u%level)

        call INS_Integrator_BDF_PrepStep_3D( ml_ins % ins_op(l), t, dt        &
                                           , u    % level(l) % val            &
                                           , u1   % level(l) % val            &
                                           , f_c1 % level(l) % val            &
                                           , f_d1 % level(l) % val            &
                                           , tau                              &
                                           , f    % level(l) % val            &
                                           , f_d  % level(l) % val            &
                                           , bv   % level(l) % var            &
                                           , mu   % level(l) % val(:,:,:,:,1) &
                                           , nu   % level(l) % val(:,:,:,:,1) &
                                           , first                            )

      end do

      ! extrapolation ..........................................................

      do l = 1, size(u%level)
        associate( mesh_l => ml_ins % ml_op_u % sem(l) % mesh  &
                 , v_l    => u   % level(l) % val(:,:,:,:,1:3) &
                 , f_l    => f   % level(l) % val(:,:,:,:,1:3) &
                 , f_d0_l => f_d % level(l) % val(:,:,:,:,1:3) )

          ! update frozen elements
          if (l > 1) then
            associate( mesh_p => ml_ins % ml_op_u % sem(l-1) % mesh &
                     , iop_pl => ml_ins % ml_op_u % iop_cf_x(l-1)   &
                     , v_p    => u % level(l-1) % val(:,:,:,:,1:3)  )

              ! restrict to frozen elements, once possible
              call ParentToChildInterpolation_3D &
                       (mesh_p, mesh_l, iop_pl, v_p, v_l)

            end associate
          end if

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, mesh_l % n_elem_active
            v_l(:,:,:,e,c) = tau * (f_l(:,:,:,e,c) + f_d0_l(:,:,:,e,c))
          end do
          end do

        end associate
      end do

      if (this % coupled) then

        ! coupled Stokes solver ................................................

        call ml_ins % Stokes_MG_Start(tau, mu, nu, bv, f_d, f, u, this%n_fmg)
        call ml_ins % Stokes_MG_Cycle(tau, mu, nu, bv, f, u, this%n_cyc)

      else

        ! projection ...........................................................

        call ml_ins % ProjectionStep(tau, bv, u)

        ! diffusion RHS ........................................................

        do l = 1, size(u%level)
          associate( mesh_l => ml_ins % ins_op(l) % mesh         &
                   , v_l    => u   % level(l) % val(:,:,:,:,1:3) &
                   , f_l    => f   % level(l) % val(:,:,:,:,1:3) &
                   , f_d0_l => f_d % level(l) % val(:,:,:,:,1:3) )

            !$omp do collapse(2)
            do c = 1, 3
            do e = 1, mesh_l % n_elem_active
              f_l(:,:,:,e,c) = 1/tau * v_l(:,:,:,e,c) - f_d0_l(:,:,:,e,c)
            end do
            end do
          end associate

        end do

        ! update outflow conditions ............................................

        block
          type(BoundaryVariable_3D), save :: bv_dp, bv_p
          integer :: b

          do b = 1, size(problem % bc_v)
            if (problem % bc_v(b) /= 'O') cycle

            do l = 1, size(u%level)

              !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call bv % level(l) % var(b) % GetSlice(bv_dp, first=3, last=3)
              call bv % level(l) % var(b) % GetSlice(bv_p , first=4, last=4)
              !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !$omp barrier

              ! pᵇ = p
              call bv_p % Extract(u % level(l) % val(:,:,:,:,4))

              ! pᵇ = p - ∆pᵇ
              call MergeArrays(  ONE, bv_p  % val(:,:,:,1) &
                              , -ONE, bv_dp % val(:,:,:,1) )

            end do
          end do
        end block

        ! viscous diffusion ....................................................

        call ml_ins % DiffusionStep(tau, mu, nu, bv, f, u)

      end if

      !$omp master
      t = t + dt
      !$omp end master

      ! cleanup ................................................................

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (last) then
        deallocate(f, f_d, mu, nu, u1, f_c1, f_d1, bv)
      end if
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'exit'
    end if

  end subroutine TimeStep

  !=============================================================================
  ! TBP of ML_INS_Integrator_BDF_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of BDF time-integrator options

  subroutine Bcast_ML_INS_Integrator_BDF_Options(this, root, comm)
    class(ML_INS_Integrator_BDF_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call this % ML_INS_IntegratorOptions_3D % Bcast(root, comm)

    call XMPI_Bcast(this % coupled, root, comm)
    call XMPI_Bcast(this % n_fmg  , root, comm)
    call XMPI_Bcast(this % n_cyc  , root, comm)

  end subroutine Bcast_ML_INS_Integrator_BDF_Options

  !=============================================================================

end module ML__INS__Integrator__BDF__3D
