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

!> summary:  2D Cylinder benchmark
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/10/11, revised 2026/11/02
!>
!> References
!>
!> 1. Schäfer M, Turek S. The benchmark problem 'flow around a cylinder’.
!>    In Flow Simulation with High-Performance Computers II, Hirschel EH (ed.).
!>    Notes Num Fluid Mech 52. Vieweg: Braunschweig 547–566, 1996
!>
!> 2. John V. Reference values for drag and lift of a two‐dimensional
!>    time‐dependent flow around a cylinder.
!>    Int J Num Meth Fluids, 44(7):777-788, 2004
!>
!> Problem 1 is the case 2D-2 with stationary inflow at Re = 100, and
!> problem 2 is the case 2D-3 with oscillating inflow at Re = 100.
!>
!===============================================================================

module INS__Problem__Cylinder_2D__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Cylinder2D_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_Cylinder2D_3D

    integer :: test_case !< 1 for steady, and 2 for unsteady flow

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type INS_Problem_Cylinder2D_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Cylinder2D_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    integer :: test_case = 1
    character, allocatable :: bc_v(:)

    namelist /parameters/ test_case, bc_v

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

    allocate(bc_v(nb), source = 'D')

    bc_v(2) = 'O'

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(test_case, 0, comm)
      call XMPI_Bcast(bc_v     , 0, comm)
    end if

    problem % test_case       =  test_case
    problem % stokes          =  .false.
    problem % exact_solution  =  .false.
    problem % v_ref           =  1.5E+0_RNP
    problem % nu_ref          =  1.0E-3_RNP

    call move_alloc(bc_v, problem % bc_v)
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Cylinder2D_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    call SetArray(u, ZERO, multi = .true.)

    if (problem % nc == 4 .or. size(x) == 0) return

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Cylinder2D_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    real(RNP), parameter :: hc = 0.41_RNP

    real(RNP) :: c0, y
    integer   :: np, nf, nc
    integer   :: i, j, k

    np = size(ub, 1)
    nf = size(ub, 3)
    nc = size(ub, 4)

    call SetArray(ub, ZERO, multi = .true.)

    if (b /= 1) return

    select case(problem % test_case)
    case(1)
      ! steady, Re = 100
      c0 = 4 * problem % v_ref / hc ** 2
    case(2)
      ! unsteady steady, Re = 100
      c0 = real(4 * problem % v_ref * sin(PI * t/8) / hc ** 2, RNP)
    end select

    !$omp do collapse(3)
    do k = 1, nf
    do j = 1, np
    do i = 1, np
      y = xb(i,j,k,2)
      ub(i,j,k,1) = c0 * y * (hc - y)
    end do
    end do
    end do

  end subroutine GetBoundaryValues

  !=============================================================================

end module INS__Problem__Cylinder_2D__3D
