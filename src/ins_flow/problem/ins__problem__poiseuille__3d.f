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

!> summary:  Poiseuille flow
!> author:   Joerg Stiller
!> date:     2020/08/11
!>
!> For Poiseuille flow the Navier-Stokes equations simplify to
!>
!>       ρν d²u/dy² = -f
!>
!> with u = v₁ and BC
!>
!>       u(0) = u(h) = 0
!>
!> where u = v₁ and y = x₂. Here we use the nondimensional form
!>
!>       1/Re d²u/dy² = -f
!>
!> based on the bulk Reynolds number `Re = u_m h / 2ν`. For simplicity we
!> further impose `ρ = 1`, `δ = h/2 = 1` and `u_m = 1` such that `ν = 1/Re`.
!> Using these assumptions, the relation for the bulk velocity
!>
!>       u_m = f h² / (12ρν)
!>
!> yields `f = 3 / Re` and the exact solution takes the form
!>
!>       u(y) = 6 (1 - y/2) y/2
!>
!===============================================================================

module INS__Problem__Poiseuille__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Poiseuille_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_Poiseuille_3D

    real(RNP) :: re     !< bulk Reynolds number, Re = 1/ν
    real(RNP) :: alpha  !< max relative perturbation of velocity components

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources

    procedure :: GetExactSolution
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_Poiseuille_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Poiseuille_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .false. ! Navier-Stokes is default
    real(RNP) :: re     = 1       ! bulk flow Reynolds number, Re = 1/ν
    real(RNP) :: alpha  = 0       ! max relative perturbation
    real(RNP) :: l      = 2       ! channel length
    real(RNP) :: w      = 2       ! channel width
    character, allocatable :: bc_v(:)

    namelist /parameters/ stokes, re, alpha, l, w, bc_v

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
      call XMPI_Bcast(stokes, 0, comm)
      call XMPI_Bcast(re    , 0, comm)
      call XMPI_Bcast(alpha , 0, comm)
      call XMPI_Bcast(l     , 0, comm)
      call XMPI_Bcast(w     , 0, comm)
      call XMPI_Bcast(bc_v  , 0, comm)
    end if

    problem % stokes          =  stokes
    problem % exact_solution  =  .true.
    problem % re              =  re
    problem % alpha           =  alpha
    problem % v_ref           =  ONE
    problem % nu_ref          =  ONE / re
    problem % x0              =  ZERO
    problem % x1              =  [l, TWO, w]

    call move_alloc(bc_v, problem % bc_v)
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Poiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,2))

    call GetVelocity( n, alpha = problem % alpha  &
                    , y = x(:,:,:,:, 2)           &
                    , v = u(:,:,:,:,1:3)          )

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Poiseuille_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,2))

    call GetVelocity(n, alpha = ZERO, y = xb(:,:,:,2), v = ub(:,:,:,1:3))

    ! remaining variables get zero
    do m = 4, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0 .or. t < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_Poiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    ! f₁ = 3 / Re
    call SetArray(F_s(:,:,:,:,1 ), 3 / problem % re)
    call SetArray(F_s(:,:,:,:,2:), ZERO, multi=.true.)

    ! silence the compiler ;)
    if (size(x) < 0 .or. t < 0) return

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t), or zero if not available

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_Poiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,2))

    call GetVelocity(n, alpha = ZERO, y = x(:,:,:,:,2), v = u(:,:,:,:,1:3))

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Provides the exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_Poiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    call SetArray(F_d(:,:,:,:,1 ), -3 / problem % re)
    call SetArray(F_d(:,:,:,:,2:), ZERO, multi=.true.)

    ! silence the compiler ;)
    if (size(x) < 0 .or. t < 0) return

  end subroutine GetExactDiffusiveTerm

  !=============================================================================
  ! Problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, alpha, y, v)
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: alpha  !< max relative perturbation
    real(RNP), intent(in)  :: y(n)   !< y-coordinate of points
    real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

    real(RNP) :: a, yh
    integer   :: i

    call SetArray(v, ZERO, multi = .true.)

    if (alpha > ZERO) then
      call random_number(v)  ! random numbers in the range [0,1)
    end if

    ! exact velocity plus random values scaled to ±α u_c
    !$omp do
    do i = 1, n

      yh = HALF * y(i)

      ! eliminate fluctuations on the wall
      if (min(yh, 1-yh) < epsilon(yh)) then
        a = 0
      else
        a = alpha
      end if

      v(i,1) = a * (2*v(i,1) - 1)  +  6 * (1 - yh) * yh
      v(i,2) = a * (2*v(i,2) - 1)
      v(i,3) = a * (2*v(i,3) - 1)

    end do

  end subroutine GetVelocity

  !=============================================================================

end module INS__Problem__Poiseuille__3D
