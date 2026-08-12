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

!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/03/05
!===============================================================================

module CART__ISP_Flow__Time_Integrator
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use XMPI
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  implicit none
  private

  public :: TimeIntegrator
  public :: TimeIntegratorOptions

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for incompressible flow
  !>
  !>   * `project`                                                           \n
  !>        0: no additional projection step                                 \n
  !>        1: additional projection at the end of the time step             \n
  !>        2: additional projection at the end of each stage (RK only)
  !>
  !>   * `pressure`
  !>        0: return pressure as computed                                   \n
  !>        1: recompute pressure at the end of the time step

  type, abstract :: TimeIntegrator

    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators

    character(len=80) :: name = ''  !< time-integrator name

    ! options
    integer :: project    = 0 !< switch for additional projections
    integer :: pressure   = 0 !< switch for pressure recomputation

  contains
    procedure, non_overridable :: Init_TimeIntegrator
    procedure, non_overridable :: Show_TimeIntegrator
    procedure :: Show => Show_TimeIntegrator
    procedure(TimeStep), deferred :: TimeStep
  end type TimeIntegrator

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step
    !>
    !> Activate the `standby` option to allow for reusing workspace and data

    subroutine TimeStep(this, t, dt, u, standby)
      import
      class(TimeIntegrator), intent(inout) :: this
      real(RNP),         intent(inout) :: t            !< time t₀ → t
      real(RNP),         intent(in)    :: dt           !< step size ∆t = t-t₀
      real(RNP),         intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
      logical, optional, intent(in)    :: standby      !< reuse workspace
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type TimeIntegratorOptions
    integer :: project    = 0 !< switch for additional projections
    integer :: pressure   = 0 !< switch for pressure recomputation
  contains
    procedure :: Bcast => Bcast_TimeIntegratorOptions
  end type TimeIntegratorOptions

contains

  !=============================================================================
  ! TimeIntegrator: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of TimeIntegrator object

  subroutine Init_TimeIntegrator(this, problem, flow_op, opt)
    class(TimeIntegrator),        intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem    !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op    !< flow operators
    class(TimeIntegratorOptions), optional, intent(in) :: opt !< options

    this % problem => problem
    this % flow_op => flow_op

    ! options
    if (present(opt)) then
      this % project    = opt % project
      this % pressure   = opt % pressure
    end if

  end subroutine Init_TimeIntegrator

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator settings

  subroutine Show_TimeIntegrator(this, unit)
    class(TimeIntegrator), intent(in) :: this
    integer,     optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)')       'TimeIntegrator settings'
    write(io,'(A,/)')       repeat('=',80)
    write(io,'(2X,A,T15,A)')  'name:'     , trim(this % name)
    write(io,'(2X,A,T15,I0)') 'project:'  , this % project
    write(io,'(2X,A,T15,I0)') 'pressure:' , this % pressure
    write(io,*)

    ! append further settings in corresponding routines of derived types

  end subroutine Show_TimeIntegrator

  !=============================================================================
  ! TimeIntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_TimeIntegratorOptions(this, root, comm)
    class(TimeIntegratorOptions), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(2)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % project  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % pressure , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_TimeIntegratorOptions

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator
