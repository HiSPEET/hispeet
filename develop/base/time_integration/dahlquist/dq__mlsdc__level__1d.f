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

!> summary:  MLSDC level for Dahlquist
!> author:   Erik Pfister
!> date:     2024/04/15
!===============================================================================

module DQ__MLSDC__Level__1D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use Execution_Control
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
  use DQ__Time_Integrator
  use DQ__SDC__Method
  use DQ__MLSDC__Corrector__Euler
  use DQ__MLSDC__Corrector__ISD

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Type for keeping one level of the MLSDC data structure

  type, public :: DQ_MLSDC_Level_1D

    logical :: is_root  !< T if root (bottom) level mesh
    logical :: is_top   !< T if top level mesh

    integer :: p_time   !< polynomial degree in time
    integer :: m_time   !< number of subintervals per time step
    integer :: n_time   !< number of steps in time
    integer :: n_stage_p
    integer :: n_stage_c

    ! discretization and solvers ...............................................

    class(DQ_SDC_Method), allocatable :: sdc

    ! transfer operators .......................................................

    type(HP_RefinementOperator_1D) :: iop_cf_t !< C-F interpolation in t

    type(HP_CoarseningOperator_1D) :: pop_fc_t !< F-C projection in t

  contains

    procedure :: GetTimeMesh    !< get the time mesh points
    procedure :: GetResidual    !< compute multi-step collocation residual
    procedure :: ApplyOperator  !< application of operator
    procedure :: ApplyPredictor !< application of SDC predictor
    procedure :: ApplyCorrector !< application of SDC corrector
    procedure :: Interpolate_CF !< solution interpolation to next finer   level
    procedure :: Project_FC     !< solution interpolation to next coarser level
    procedure :: Restrict_FC    !< residual restriction from next finer   level

  end type DQ_MLSDC_Level_1D

  ! constructor
  interface DQ_MLSDC_Level_1D
    module procedure New_DQ_MLSDC_Level_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for DQ_MLSDC_Level_1D initialization

  type, public :: DQ_MLSDC_Level_Options_1D

    integer :: p_time  = -1 !< polynomial degree in time
    integer :: n_time  = -1 !< number of steps in time
    integer :: n_stage_p = -1
    integer :: n_stage_c = -1

    class(DQ_TimeIntegrator_Options), allocatable :: pre_opt
    class(DQ_SDC_Options),            allocatable :: sdc_opt

    type(HP_RefinementOptions_1D) :: iop_cf_t

    type(HP_CoarseningOptions_1D) :: pop_fc_t

  end type DQ_MLSDC_Level_Options_1D

  !=============================================================================
  ! module procedures

  interface

    !---------------------------------------------------------------------------
    !> Application of the SDC predictor

    module subroutine ApplyPredictor(this, lambda, dt, u)
      class(DQ_MLSDC_Level_1D), intent(inout) :: this
      complex(RNP), intent(in)    :: lambda  !< lambda
      real(RNP)   , intent(in)    :: dt      !< size of the time slice
      complex(RNP), intent(inout) :: u(0:,:) !< approximate solution
    end subroutine ApplyPredictor

    !---------------------------------------------------------------------------
    !> Application of the SDC corrector

    module subroutine ApplyCorrector( &
        this, incremental, lambda, dt, G, u, n_sweep)
      class(DQ_MLSDC_Level_1D), intent(inout) :: this
      logical     , intent(in)    :: incremental
      complex(RNP), intent(in)    :: lambda  !< lambda
      real(RNP)   , intent(in)    :: dt      !< size of the time slice
      complex(RNP), intent(in)    :: G(0:,:) !< FAS defect correction
      complex(RNP), intent(inout) :: u(0:,:) !< approximate solution
      integer     , intent(in)    :: n_sweep !< number of sweeps
    end subroutine ApplyCorrector

    !---------------------------------------------------------------------------
    !> Application of the operator

    module subroutine ApplyOperator(this, incremental, lambda, dt, u, r)
      class(DQ_MLSDC_Level_1D), intent(in) :: this
      logical     , intent(in)  :: incremental
      complex(RNP), intent(in)  :: lambda  !< lambda
      real(RNP),    intent(in)  :: dt      !< size of the time slice
      complex(RNP), intent(in)  :: u(0:,:) !< approximate solution
      complex(RNP), intent(out) :: r(0:,:) !< u after operator applied
    end subroutine ApplyOperator

    !---------------------------------------------------------------------------
    !> Residual of the collocation method

    module subroutine GetResidual(this, incremental, lambda, dt, G, u, r)
      class(DQ_MLSDC_Level_1D), intent(in) :: this !< MLSDC level
      logical,      intent(in)  :: incremental
      complex(RNP), intent(in)  :: lambda  !< lambda
      real(RNP),    intent(in)  :: dt      !< size of the time slice
      complex(RNP), intent(in)  :: G(0:,:) !< FAS RHS
      complex(RNP), intent(in)  :: u(0:,:) !< approximate solution
      complex(RNP), intent(out) :: r(0:,:) !< residual
    end subroutine GetResidual

    !---------------------------------------------------------------------------
    !> Coarse-to-fine space-time interpolation of solution-like variables

    module subroutine Interpolate_CF(this, u_c, u_f, complete)
      class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
      complex(RNP), intent(in)    :: u_c(0:,:) !< coarse solution variable
      complex(RNP), intent(inout) :: u_f(0:,:) !< fine solution variable
      logical,      intent(in)    :: complete  !< T/F for all/refined elements
    end subroutine Interpolate_CF

    !---------------------------------------------------------------------------
    !> Fine-to-coarse space-time projection of solution-like variables

    module subroutine Project_FC(this, u_f, u_c)
      class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
      complex(RNP), intent(in)    :: u_f(0:,:) !< fine solution variable
      complex(RNP), intent(inout) :: u_c(0:,:) !< coarse solution variable
    end subroutine Project_FC

    !---------------------------------------------------------------------------
    !> Fine-to-coarse restriction of residual-like variables

    module subroutine Restrict_FC(this, r_f, r_c)
      class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
      complex(RNP), intent(in)    :: r_f(0:,:) !< fine residual variable
      complex(RNP), intent(inout) :: r_c(0:,:) !< coarse residual variable
    end subroutine Restrict_FC

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Returns a new DQ_MLSDC_Level_1D object

  function New_DQ_MLSDC_Level_1D(opt) result(this)
    class(DQ_MLSDC_Level_Options_1D), intent(in) :: opt
    type(DQ_MLSDC_Level_1D) :: this

    call Init_DQ_MLSDC_Level_1D(this, opt)

  end function New_DQ_MLSDC_Level_1D

  !-----------------------------------------------------------------------------
  !> Inititialization of an MLSDC level

  subroutine Init_DQ_MLSDC_Level_1D(this, opt)
    class(DQ_MLSDC_Level_1D),         intent(inout) :: this
    class(DQ_MLSDC_Level_Options_1D), intent(in)    :: opt

    ! safeguard ................................................................

    if (.not. allocated(opt % pre_opt)) then
      call Error( 'Init_DQ_MLSDC_Level_1D', 'opt % pre_opt not allocated')
    end if

    if (.not. allocated(opt % sdc_opt)) then
      call Error( 'Init_DQ_MLSDC_Level_1D', 'opt % sdc_opt not allocated')
    end if

    ! parameters ...............................................................

    this % is_root = opt % pop_fc_t % po_c <= 0
    this % is_top  = opt % iop_cf_t % po_f <= 0

    this % p_time  = opt % p_time
    this % n_time  = opt % n_time
    this % n_stage_p = opt % n_stage_p
    this % n_stage_c = opt % n_stage_c

    ! discretization and solvers ...............................................

    select type(sdc_opt => opt % sdc_opt)
    type is(DQ_MLSDC_Corrector_Options_Euler)
      this % sdc = DQ_MLSDC_Corrector_Euler(opt % pre_opt, sdc_opt)
    type is(DQ_MLSDC_Corrector_Options_ISD)
      this % sdc = DQ_MLSDC_Corrector_ISD(opt % pre_opt, sdc_opt)
    end select

    this % m_time = this % sdc % n_sub

    ! transfer operators .......................................................

    if (.not. this % is_top) then
      this % iop_cf_t = HP_RefinementOperator_1D(opt % iop_cf_t)
    end if

    if (.not. this % is_root) then
      this % pop_fc_t = HP_CoarseningOperator_1D(opt % pop_fc_t)
    end if

  end subroutine Init_DQ_MLSDC_Level_1D

  !-----------------------------------------------------------------------------
  !> Get the time mesh points

  subroutine GetTimeMesh(this, t_0, t_1, t)
    class(DQ_MLSDC_Level_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: t_0    !< start time
    real(RNP),                intent(in)  :: t_1    !< end time
    real(RNP), allocatable,   intent(out) :: t(:,:) !< time mesh points

    real(RNP) :: dt
    integer   :: n

    associate( m_time => this % m_time &
             , n_time => this % n_time &
             , sdc => this % sdc )

      allocate(t(0:m_time,1:n_time))

      dt = (t_1 - t_0) / n_time

      do n = 1, n_time
        t(0:,n) = sdc % SubintervalPoints(t_0 + (n-1)*dt, dt)
      end do

    end associate

  end subroutine GetTimeMesh

  !=============================================================================

end module DQ__MLSDC__Level__1D
