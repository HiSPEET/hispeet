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

!> summary:  Placeholder for a user-defined incompressible flow problem
!> author:   Joerg Stiller
!> date:     2026/10/29
!>
!> This dummy is compiled when no user-defined problem is selected. It
!> provides the module and type expected by Set_INS_TestProblem_3D for the
!> problem name 'User' and stops with an error when the problem is set.
!>
!> A user-defined problem replaces this file at configuration time:
!>
!>     cmake -D INS_USER_PROBLEM=<name>|<directory> ...
!>
!> where <name> refers to user/ins_problem/<name>. See
!> user/ins_problem/README.md for details.
!===============================================================================

module INS__Problem__User__3D

  use Kind_Parameters, only: RNP
  use Execution_Control
  use MPI_Binding

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_User_3D

  !-----------------------------------------------------------------------------
  !> Dummy type for a user-defined problem

  type, extends(INS_Problem_3D) :: INS_Problem_User_3D
  contains
    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
  end type INS_Problem_User_3D

  character(len=*), parameter :: MODULE_NAME = 'INS__Problem__User__3D'

contains

  !-----------------------------------------------------------------------------
  !> Initialization: stops with an error, since no user problem is available

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_User_3D), intent(inout) :: problem
    integer,                    intent(in)    :: nb    !< number of boundaries
    character(len=*), optional, intent(in)    :: file  !< input file
    type(MPI_Comm),   optional, intent(in)    :: comm  !< MPI communicator

    call Error( 'SetProblem'                                        &
              , 'no user-defined problem, configure HiSPEET with '  &
                // '-D INS_USER_PROBLEM=<name>|<directory>'         &
              , MODULE_NAME                                         )

    ! silence the compiler ;)
    if (nb < 0 .or. present(file) .or. present(comm)) return

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Initial values: not available

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_User_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    call Error('GetInitialValues', 'no user-defined problem', MODULE_NAME)

    u = 0

    ! silence the compiler ;)
    if (size(x) < 0) return

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Boundary values: not available

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_User_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    call Error('GetBoundaryValues', 'no user-defined problem', MODULE_NAME)

    ub = 0

    ! silence the compiler ;)
    if (b < 0 .or. size(xb) < 0 .or. t < 0) return

  end subroutine GetBoundaryValues

  !=============================================================================

end module INS__Problem__User__3D
