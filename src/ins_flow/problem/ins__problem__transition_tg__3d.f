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

!> summary:  Nonlinear transition and decay of a disturbed Taylor-Green vortex
!> author:   Joerg Stiller
!> date:     2020/08/11
!>
!> References
!>
!>   -  ME Brachet, Direct simulation of three-dimensional turbulence in the
!>      Taylor–Green vortex. Fluid Dyn Res 8(1-4):1-8, 1991
!>
!>   -  M Guesmi, M Grotteschi, J Stiller, Assessment of high-order IMEX methods
!>      for incompressible flow. Int J Numer Meth Fluids 95:954-978, 2023
!>
!===============================================================================

module INS__Problem__Transition_TG__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_TransitionTG_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_TransitionTG_3D
    real(RNP) :: nu !< kinematic viscosity
  contains
    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources
  end type INS_Problem_TransitionTG_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_TransitionTG_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    real(RNP) :: nu = 0.01 ! kinematic viscosity, ν = 1/Re
    character, allocatable :: bc_v(:)

    namelist /parameters/ nu

    logical :: exists
    integer :: prm, rank

    ! preliminaries ............................................................

    if (present(comm)) then
      call MPI_Comm_rank(comm, rank)
    else
      rank = 0
    end if

    ! check for input file
    if (rank == 0 .and. present(file)) then

      exists = len_trim(file) > 0
      if (exists) then
        inquire(file=trim(file)//'.prm', exist=exists)
      end if

      if (exists) then
        open(newunit=prm, file=trim(file)//'.prm')
      else
        call Warning('SetProblem','Input file "'//trim(file)//'.prm" not found')
      end if

    else
      exists = .false.
    end if

    ! parameters ...............................................................

    ! default BC (pressure is ignored)
    allocate(bc_v(nb), source = 'D')

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(nu  , 0, comm)
      call XMPI_Bcast(bc_v, 0, comm)
    end if

    problem % stokes         = .false.
    problem % exact_solution = .false.

    problem % v_ref  =  2   ! initial maximum velocity
    problem % nu_ref =  nu
    problem % x0     = -PI
    problem % x1     =  PI

    call move_alloc(bc_v, problem % bc_v)
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_TransitionTG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetInitialVelocity(problem, n, x, v = u(:,:,:,:,1:3))

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_TransitionTG_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    call Warning( 'GetBoundaryValues'               &
                , 'No boundary values available'    &
                , 'ISP_Flow_Problem__Transition_TG' )

    call SetArray(ub, ZERO, multi=.true.)

    ! silence the compiler ;)
    if (b < 0 .or. size(xb) < 0 .or. t < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_TransitionTG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    call SetArray(F_s, ZERO, multi=.true.)

    ! silence the compiler ;)
    if (size(x) < 0 .or. t < 0) return

  end subroutine GetExternalSources

  !=============================================================================
  ! Problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Initial velocity

  subroutine GetInitialVelocity(problem, n, x, v)
    class(INS_Problem_TransitionTG_3D), intent(in) :: problem
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: x(n,3) !< mesh points
    real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

    real(RNP) :: x1, x2, x3 ! coordinates xi
    integer   :: i

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      x3 = x(i,3)
      v(i,1) =  cos(x3) * sin(x1) * cos(x2)
      v(i,2) = -cos(x3) * sin(x2) * cos(x1)
      v(i,3) =  ZERO
    end do

  end subroutine GetInitialVelocity

  !=============================================================================

end module INS__Problem__Transition_TG__3D
