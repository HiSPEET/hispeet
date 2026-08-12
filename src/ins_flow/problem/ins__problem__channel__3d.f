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

!> summary:  3D turbulent channel flow
!> author:   Joerg Stiller
!> date:     2020/08/11
!>
!> Problem definition for turbulent channel flow
!>
!>   - flow in x-direction
!>   - x: length l, periodic
!>   - y: width w, periodic
!>   - z: walls at z = 0 and z = 2δ
!>   - volume force in x-direction set to match prescribed Re_τ = u_τ δ / ν
!>   - initial conditions (may be improved ;)
!>       + parabolic profile matching expected mean velocity ū
!>       + random fluctuation with amplitude αū
!>
!> References
!>
!>   -  J Kim, P Moin, RD Moser, Turbulence statistics in fully developed
!>      turbulent channel flow at low Reynolds number. JFM 177:133-166, 1987
!>
!>   -  KS Nelson, OB Fringer, Reducing spin-up time for simulations of
!>      turbulent channel flow. Phys Fluids 29:105101, 2017
!>
!===============================================================================

module INS__Problem__Channel__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Channel_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_Channel_3D

    real(RNP) :: re_t   !< friction Reynolds number, Re_τ = u_τ δ / ν
    real(RNP) :: delta  !< channel half height δ
    real(RNP) :: alpha  !< max relative perturbation of velocity components

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources

  end type INS_Problem_Channel_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Channel_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    ! problem parameters according to Kim et al. (1987)
    real(RNP) :: re_t  = 180  ! friction Reynolds number, Re_τ = u_τ δ / ν
    real(RNP) :: alpha = 1    ! max relative perturbation of velocity components
    real(RNP) :: l     = 4*PI ! channel length
    real(RNP) :: h     = 2    ! channel height, h = 2δ
    real(RNP) :: w     = 2*PI ! channel width
    character, allocatable :: bc_v(:)

    namelist /parameters/ re_t, alpha, l, h, w

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
      call XMPI_Bcast(re_t , 0, comm)
      call XMPI_Bcast(alpha, 0, comm)
      call XMPI_Bcast(l    , 0, comm)
      call XMPI_Bcast(h    , 0, comm)
      call XMPI_Bcast(w    , 0, comm)
      call XMPI_Bcast(bc_v , 0, comm)
    end if

    problem % stokes          =  .false.
    problem % exact_solution  =  .false.
    problem % re_t            =  re_t
    problem % delta           =  h / 2
    problem % alpha           =  alpha
    problem % v_ref           =  14.64 / h * re_t**(ONE/7) ! v_bulk (Dean 1978)
    problem % nu_ref          =  ONE / re_t
    problem % x0              =  ZERO
    problem % x1              = [l, w, h]

    call move_alloc(bc_v, problem % bc_v)
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Channel_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetInitialVelocity(problem, n, z = x(:,:,:,:,3), v = u(:,:,:,:,1:3))

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Channel_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    call SetArray(ub, ZERO, multi=.true.)

    ! silence the compiler ;)
    if (b < 0 .or. size(xb) < 0 .or. t < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_Channel_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    ! f₁ = (u_τ)² / δ
    call SetArray(F_s(:,:,:,:,1 ), ONE / problem % delta ** 3)
    call SetArray(F_s(:,:,:,:,2:), ZERO, multi=.true.)

    ! silence the compiler ;)
    if (size(x) < 0 .or. t < 0) return

  end subroutine GetExternalSources

  !=============================================================================
  ! Problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Initial velocity

  subroutine GetInitialVelocity(problem, n, z, v)
    class(INS_Problem_Channel_3D), intent(in) :: problem
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: z(n)   !< z-coordinate of points
    real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

    real(RNP) :: a, zh
    integer   :: i

    associate( u_m  => problem % v_ref )

      call SetArray(v, ZERO, multi = .true.)

      ! initialize v with random numbers in the range [0,1)
      call random_seed()     ! initialize with system generated seed
      call random_number(v)

      ! rescale random values to ±α u_m and add quadratic profile to v₁
      !$omp do
      do i = 1, n

        zh = z(i) / (2 * problem % delta)

        ! eliminate fluctuations on the wall
        if (min(zh, 1-zh) < epsilon(zh)) then
          a = 0
        else
          a = problem % alpha * u_m
        end if

        v(i,1) = a * (2*v(i,1) - 1)  +  6 * u_m * (1 - zh) * zh
        v(i,2) = a * (2*v(i,2) - 1)
        v(i,3) = a * (2*v(i,3) - 1)

      end do

    end associate

  end subroutine GetInitialVelocity

  !=============================================================================

end module INS__Problem__Channel__3D
