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

!> summary:  No-flow test case extended from Akbas et al. (2018)
!> author:   Joerg Stiller
!> date:     2024/03/22
!>
!> Problem intended to show the necessity of grad-div stabilization.
!> The pressure is prescribed as follows
!>   - in the 2D case, Ω = (0,1)²:  p = sin 2π(x+y)
!>   - in the 3D case, Ω = (0,1)³:  p = sin 2π(x+y+y)
!>
!> and the force distribution f = ∇p. The velocity is v = 0.
!>
!> ### References
!>
!> 1. Akbas M, Linke A, Rebholz LG & Schroeder PW, Comput Meth Appl Mech Engrg
!>    341:917-938, 2018
!===============================================================================

module INS__Problem__No_Flow__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_NoFlow_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the No-Flow problem

  type, extends(INS_Problem_3D) :: INS_Problem_NoFlow_3D
    integer :: dim !< dimension of the problem {2,3}
  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources

    procedure :: GetExactSolution
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_NoFlow_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_NoFlow_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .true. ! set true for Stokes problem
    integer   :: dim    = 2      ! select 2D or 3D problem {2,3}
    real(RNP) :: nu     = 1e-4   ! kinematic viscosity
    real(RNP) :: x0(3)  = 0      ! bounding box corner nearest to to -∞
    real(RNP) :: x1(3)  = 1      ! bounding box corner nearest to to +∞

    namelist /parameters/ stokes, dim, nu, x0, x1

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

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(stokes, 0, comm)
      call XMPI_Bcast(dim   , 0, comm)
      call XMPI_Bcast(nu    , 0, comm)
      call XMPI_Bcast(x0    , 0, comm)
      call XMPI_Bcast(x1    , 0, comm)
    end if

    problem % stokes         = stokes
    problem % exact_solution = .true.
    problem % dim            = dim
    problem % nu_ref         = nu
    problem % v_ref          = ONE
    problem % x0             = x0
    problem % x1             = x1

    allocate(problem % bc_v(nb), source = 'D')
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    call SetArray(u, ZERO, multi = .true.)
    call GetPressure(problem%dim, size(x(:,:,:,:,1)), x, u(:,:,:,:,4))

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` in points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    call SetArray(ub, ZERO, multi=.true.)
    call GetPressure(problem%dim, size(xb(:,:,:,1)), xb, ub(:,:,:,4))

    ! ignore boundary identifier and time
    if (b < 0 .or. t > 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    call GetPressureTerm(problem%dim, size(x(:,:,:,:,1)), x, F_s)
    call ScaleArray(F_s, -ONE, multi=.true.)

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    call SetArray(u, ZERO, multi = .true.)
    call GetPressure(problem%dim, size(x(:,:,:,:,1)), x, u(:,:,:,:,4)  )

    if (t > 0) return

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< pressure term

    call GetPressureTerm(problem%dim, size(x(:,:,:,:,1)), x, F_p)

    if (t > 0) return

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_NoFlow_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    call SetArray(F_d, ZERO, multi=.true.)

    if (t > 0) return

  end subroutine GetExactDiffusiveTerm

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(dim, n, x, p)
    integer,   intent(in)  :: dim     !< problem dimension
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    integer :: i

    select case(dim)
    case(2)
      !$omp do
      do i = 1, n
        p(i) = sin(2*PI * (x(i,1) + x(i,2)))
      end do
    case default
      !$omp do
      do i = 1, n
        p(i) = sin(2*PI * (x(i,1) + x(i,2) + x(i,3)))
      end do
    end select

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure term

  subroutine GetPressureTerm(dim, n, x, F_p)
    integer,   intent(in)  :: dim       !< problem dimension
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(out) :: F_p(n,3)  !< pressure term, -∇p

    integer :: i

    select case(dim)
    case(2)
      !$omp do
      do i = 1, n
        F_p(i,1)  = -2*PI * cos(2*PI * (x(i,1) + x(i,2)))
        F_p(i,2)  = F_p(i,1)
        F_p(i,3:) = ZERO
      end do
    case default
      !$omp do
      do i = 1, n
        F_p(i,1)  = -2*PI * cos(2*PI * (x(i,1) + x(i,2) + x(i,3)))
        F_p(i,2)  = F_p(i,1)
        F_p(i,3)  = F_p(i,1)
        F_p(i,4:) = ZERO
      end do
    end select

  end subroutine GetPressureTerm

  !=============================================================================

end module INS__Problem__No_Flow__3D
