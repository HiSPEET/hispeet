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

!> summary:  Method of Manufactured Solution for vortex field with different
!>           cases of variable diffusivity
!> author:   Karl Schoppmann, Jörg Stiller
!> date:     2019/05/02
!>
!>### Method of Manufactured Solution for vortex field with different
!>    cases of variable diffusivity
!>
!> Provides procedures for setting
!>
!>   *  the problem
!>
!> and provides procedures for getting
!>
!>   *  initial values
!>   *  boundary values
!>   *  the time derivative of boundary values
!>   *  external sources
!>   *  an exact solution
!>   *  the time derivative of the exact solution
!>   *  a variable diffusivity
!>
!> Four cases of variable viscosity are considered:
!>
!>   1. Position dependent viscosity `nu(x)`,
!>      implementation follows
!>      Niemann, M.: Bouyancy effects in turbulent liquid metal flow, 2017
!>   2. Position and time dependent viscosity `nu(x,t)`,
!>      implementation based on case 1 with further development
!>   3. Regularized Smagorinsky model based on exact solution `nu(v_ex)`,
!>   4. Nonlinear viscosity based on kinetic energy.
!>
!> For more details [see](http://arxiv.org/abs/2001.11902) and
!> Karl Schoppmann, Verifikation von Methoden hoher Ordnung für die numerische
!> Simulation inkompressibler Strömungen mit variabler Viskosität,
!> Projektarbeit zum Forschungspraktikum, ISM, TU Dresden, 2019.
!===============================================================================

module ISP_Flow_Problem__Variable_Viscosity

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_VariableViscosity

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the vortex field problem

  type, extends(FlowProblem) :: FlowProblem_VariableViscosity

    integer   :: test_case = 1    !< variable viscosity case 1|2|3|4
    real(RNP) :: nu_0      = ONE  !< constant base diffusivity
    real(RNP) :: nu_1      = ONE  !< coefficient of variable diffusivity

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetDiffusivity

  end type FlowProblem_VariableViscosity

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization
  !>
  !> The diffusivity is split into a constant base part and a variable part of
  !> diffusivity increase. For case 3 the coefficient of variable diffusivity
  !> can be interpreted as
  !>     nu_1 = ( C_S * L_filter )**2
  !> where `C_S` denotes the Smagorinsky constant and `L_filter` the LES
  !> filter width.

  subroutine SetProblem(problem, file, comm)
    class(FlowProblem_VariableViscosity), intent(inout) :: problem
    character(len=*), optional, intent(in) :: file !< input file
    type(MPI_Comm),   optional, intent(in) :: comm !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes    = .false.
    integer   :: test_case =  1         ! 1|2|3 (`nu(x)`|`nu(x,t)`|`nu(x,t,u)`)
    real(RNP) :: nu_0      =  ONE       ! constant base diffusivity
    real(RNP) :: nu_1      =  ONE       ! coefficient of variable diffusivity
    real(RNP) :: x0(3)     = -HALF      ! bounding box: corner nearest to -∞
    real(RNP) :: x1(3)     =  HALF      ! bounding box: corner nearest to +∞
    character, allocatable :: bc(:,:)

    namelist /parameters/ stokes, test_case, nu_0, nu_1, x0, x1, bc

    logical :: exists
    integer :: prm, rank
    integer :: b1, b2, d

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
    allocate(bc(6, problem % nc), source = 'P')

    if (rank == 0) then

      ! read parameters
      if (exists) then
        read(prm, nml=parameters)
        close(prm)
      end if

      ! normalize turbulent viscosity to nu_1
      select case(test_case)
      case(3)
        nu_1 = nu_1 / (48 * PI**2)
      case(4)
        nu_1 = nu_1 / 4
      end select

    end if

    if (present(comm)) then
      call XMPI_Bcast(stokes   , 0, comm)
      call XMPI_Bcast(test_case, 0, comm)
      call XMPI_Bcast(nu_0     , 0, comm)
      call XMPI_Bcast(nu_1     , 0, comm)
      call XMPI_Bcast(x0       , 0, comm)
      call XMPI_Bcast(x1       , 0, comm)
      call XMPI_Bcast(bc       , 0, comm)
    end if

    problem % stokes = stokes
    problem % exact_solution = .true.
    problem % variable_props = .true.

    allocate(problem % nu_ref( problem%nc ), source = ZERO)
    problem % nu_ref(1:3) = nu_0 + HALF * nu_1

    problem % test_case = test_case
    problem % nu_0      = nu_0
    problem % nu_1      = nu_1
    problem % x0        = x0
    problem % x1        = x1

    ! pressure BC
    do d = 1, 3
      b1 = 2*d - 1
      b2 = b1  + 1
      if (all(bc(b1:b2, d) == 'P')) then
        bc(b1:b2, 4) = 'P'
      else
        bc(b1:b2, 4) = 'N'
      end if
    end do

    call move_alloc(bc, problem % bc)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values `u(x,0)` for mesh points `x`

  subroutine GetInitialValues(problem, x, u)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, ZERO, u(:,:,:,:,1:3))
    call GetPressure(n, x, ZERO, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    call GetVelocity(n, xb, t, ub(:,:,:,1:3))
    call GetPressure(n, xb, t, ub(:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    integer,   intent(in)  :: b              !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: dt_ub(:,:,:,:) !< time derivative on boundary

    integer :: m, n

    n = size(xb(:,:,:,1))

    call GetVelocityTimeDerivative(n, xb, t, dt_ub(:,:,:,1:3))
    call GetPressureTimeDerivative(n, xb, t, dt_ub(:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(dt_ub,4)
      call SetArray(dt_ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0) return

  end subroutine GetBoundaryTimeDerivative

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points `x` and time `t`

  subroutine GetExternalSources(problem, x, t, f)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call Get_MMS_Source(problem, n, x, t, f(:,:,:,:,1:3))

    do m = 4, size(f,5)
      call SetArray(f(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExternalSources

  !-----------------------------------------------------------------------------
  !> Provides the exact solution `u(x,t)`, or zero if not available

  subroutine GetExactSolution(problem, x, t, u)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: u(:,:,:,:,:)   !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, t, u(:,:,:,:,1:3))
    call GetPressure(n, x, t, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Provides the time derivative of the exact solution, `∂u/∂t(x,t)`

  subroutine GetExactTimeDerivative(problem, x, t, dt_u)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocityTimeDerivative(n, x, t, dt_u(:,:,:,:,1:3))
    call GetPressureTimeDerivative(n, x, t, dt_u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(dt_u,5)
      call SetArray(dt_u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactTimeDerivative

  !-----------------------------------------------------------------------------
  !> Provides a variable diffusivity depending on input parameter `test_case`:
  !> 1. `nu(x)`
  !> 2. `nu(x,t)`
  !> 3. `nu(x,t,u_ex)`, where `u_ex` is the exact solution
  !> 4. `nu(x,t,u)`   , where `u`    is the approximate solution

  subroutine GetDiffusivity(problem, x, t, u, nu)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)     :: x (:,:,:,:,:) !< mesh points
    real(RNP), intent(in)     :: t             !< time
    real(RNP), intent(in)     :: u (:,:,:,:,:) !< flow variables
    real(RNP), intent(out)    :: nu(:,:,:,:,:) !< diffusivity

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu_0 => problem % nu_0, nu_1 => problem % nu_1,        &
              test_case => problem % test_case                       )
      select case(test_case)
        case(1)
          call GetDiffusivity_Case_1(n, x, nu_0, nu_1, nu)     ! nu(x)
        case(2)
          call GetDiffusivity_Case_2(n, x, t, nu_0, nu_1, nu)  ! nu(x,t)
        case(3)
          call GetDiffusivity_Case_3(n, x, t, nu_0, nu_1, nu)  ! nu(x,t,u_ex)
        case(4)
          call GetDiffusivity_Case_4(n, nu_0, nu_1, u, nu)     ! nu(x,t,u)
      end select
    end associate

    do m = 4, size(nu,5)
      call SetArray(nu(:,:,:,:,m), ZERO)
    end do

  end subroutine GetDiffusivity

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, x, t, v)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      v(i,1) = ( sin(kxt) + cos(kyt) ) * sin(kzt)
      v(i,2) = ( cos(kxt) + sin(kyt) ) * sin(kzt)
      v(i,3) = ( cos(kxt) + cos(kyt) ) * cos(kzt)
    end do

  end subroutine GetVelocity

  !-----------------------------------------------------------------------------
  !> Velocity time derivative

  subroutine GetVelocityTimeDerivative(n, x, t, dt_v)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(out) :: dt_v(n,3)  !< time derivative of velocity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)

      dt_v(i,1) =  k * ( ( cos(kxt) - sin(kyt)) * sin(kzt) &
                       + ( sin(kxt) + cos(kyt)) * cos(kzt) )
      dt_v(i,2) =  k * ( (-sin(kxt) + cos(kyt)) * sin(kzt) &
                       + ( cos(kxt) + sin(kyt)) * cos(kzt) )
      dt_v(i,3) = -k * ( ( sin(kxt) + sin(kyt)) * cos(kzt) &
                       + ( cos(kxt) + cos(kyt)) * sin(kzt) )
    end do

  end subroutine GetVelocityTimeDerivative

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(n, x, t, p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      p(i) = sin(kxt) * sin(kyt) * sin(kzt)
    end do

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure time derivative

  subroutine GetPressureTimeDerivative(n, x, t, dt_p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: dt_p(n) !< time derivative of pressure

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      dt_p(i) = k * ( cos(kxt) * sin(kyt) * sin(kzt) &
                    + sin(kxt) * cos(kyt) * sin(kzt) &
                    + sin(kxt) * sin(kyt) * cos(kzt) )
    end do

  end subroutine GetPressureTimeDerivative

  !-----------------------------------------------------------------------------
  !> Diffusivity, case 1: position dependent diffusivity `nu(x)`

  subroutine GetDiffusivity_Case_1(n, x, nu_0, nu_1, nu)
    integer,   intent(in)  :: n        !< number of points
    real(RNP), intent(in)  :: x(n,3)   !< mesh points
    real(RNP), intent(in)  :: nu_0     !< constant base diffusivity
    real(RNP), intent(in)  :: nu_1     !< coefficient of variable diffusivity
    real(RNP), intent(out) :: nu(n,3)  !< kinematic viscosity

    real(RNP) :: k, kx, ky, kz
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kx = k * x(i,1)
      ky = k * x(i,2)
      kz = k * x(i,3)

      nu(i,1) = nu_0 + nu_1 * sin(kx)**2 * sin(ky)**2 * sin(kz)**2
      nu(i,2) = nu(i,1)
      nu(i,3) = nu(i,1)
    end do

  end subroutine GetDiffusivity_Case_1

  !-----------------------------------------------------------------------------
  !> Diffusivity, case 2: position and time dependent diffusivity `nu(x,t)`

  subroutine GetDiffusivity_Case_2(n, x, t, nu_0, nu_1, nu)
    integer,   intent(in)  :: n        !< number of points
    real(RNP), intent(in)  :: x(n,3)   !< mesh points
    real(RNP), intent(in)  :: t        !< time
    real(RNP), intent(in)  :: nu_0     !< constant base diffusivity
    real(RNP), intent(in)  :: nu_1     !< coefficient of variable diffusivity
    real(RNP), intent(out) :: nu(n,3)  !< kinematic viscosity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k    = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) - t)
      kyt = k * (x(i,2) - t)
      kzt = k * (x(i,3) - t)

      nu(i,1) = nu_0 + nu_1 * sin(kxt)**2 * sin(kyt)**2 * sin(kzt)**2
      nu(i,2) = nu(i,1)
      nu(i,3) = nu(i,1)
    end do

  end subroutine GetDiffusivity_Case_2

  !-----------------------------------------------------------------------------
  !> Diffusivity, case 3: position, time and solution dependent diffusivity
  !> `nu(x,t,u)` whose variable part is based on the Smagorinsky turbulence
  !> model.

  subroutine GetDiffusivity_Case_3(n, x, t, nu_0, nu_1, nu)
    integer,   intent(in)  :: n        !< number of points
    real(RNP), intent(in)  :: x(n,3)   !< mesh points
    real(RNP), intent(in)  :: t        !< time
    real(RNP), intent(in)  :: nu_0     !< constant base diffusivity
    real(RNP), intent(in)  :: nu_1     !< coefficient of variable diffusivity
    real(RNP), intent(out) :: nu(n,3)  !< kinematic viscosity

    real(RNP), allocatable, save :: strain_rate(:,:,:)
    integer :: i, l, m

    !$omp single
    allocate(strain_rate(n,3,3))
    !$omp end single

    call GetStrainRate(n, x, t, strain_rate)

    !$omp do
    do i = 1, n
      nu(i,1) = ZERO
    end do

    do l = 1, 3
    do m = 1, 3
      !$omp do
      do i = 1, n
        ! using nu(i,1) for the square Frobenius-Norm of strain rate tensor
        nu(i,1) = nu(i,1) + strain_rate(i,m,l)**2
      end do
    end do
    end do

    !$omp do
    do i = 1, n
      ! modification of variable part (~ squared Smagorinsky-diffusivity)
      nu(i,1) = nu_0 + nu_1 * 2 * nu(i,1)
      nu(i,2) = nu(i,1)
      nu(i,3) = nu(i,1)
    end do

    !$omp barrier
    !$omp master
    deallocate(strain_rate)
    !$omp end master

  end subroutine GetDiffusivity_Case_3

  !-----------------------------------------------------------------------------
  !>

  subroutine GetDiffusivity_Case_4(n, nu_0, nu_1, v, nu)
    integer,   intent(in)  :: n        !< number of points
    real(RNP), intent(in)  :: nu_0     !< constant base diffusivity
    real(RNP), intent(in)  :: nu_1     !< coefficient of variable diffusivity
    real(RNP), intent(in)  :: v(n,3)   !< velocity
    real(RNP), intent(out) :: nu(n,3)  !< kinematic viscosity

    integer :: i

    !$omp do
    do i = 1, n
      nu(i,1) = nu_0 + nu_1 * (v(i,1)**2 + v(i,2)**2 + v(i,3)**2)
      nu(i,2) = nu(i,1)
      nu(i,3) = nu(i,1)
    end do

  end subroutine GetDiffusivity_Case_4

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient, case 1

  subroutine GetDiffusivityGradient_Case_1(n, x, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable diffusivity
    real(RNP), intent(out) :: grad_nu(n,3) !< diffusivity gradient

    real(RNP) :: k, kx, ky, kz
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kx = k * x(i,1)
      ky = k * x(i,2)
      kz = k * x(i,3)

      grad_nu(i,1) = 2 * k * nu_1 * sin(kx) * cos(kx) * sin(ky)**2 * sin(kz)**2
      grad_nu(i,2) = 2 * k * nu_1 * sin(ky) * cos(ky) * sin(kx)**2 * sin(kz)**2
      grad_nu(i,3) = 2 * k * nu_1 * sin(kz) * cos(kz) * sin(kx)**2 * sin(ky)**2
    end do

  end subroutine GetDiffusivityGradient_Case_1

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient, case 2

  subroutine GetDiffusivityGradient_Case_2(n, x, t, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable diffusivity
    real(RNP), intent(out) :: grad_nu(n,3) !< diffusivity gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k    = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) - t)
      kyt = k * (x(i,2) - t)
      kzt = k * (x(i,3) - t)

      grad_nu(i,1) = 2 * k * nu_1 * sin(kxt) * cos(kxt) &
                     * sin(kyt)**2 * sin(kzt)**2
      grad_nu(i,2) = 2 * k * nu_1 * sin(kyt) * cos(kyt) &
                     * sin(kxt)**2 * sin(kzt)**2
      grad_nu(i,3) = 2 * k * nu_1 * sin(kzt) * cos(kzt) &
                     * sin(kxt)**2 * sin(kyt)**2
    end do

  end subroutine GetDiffusivityGradient_Case_2

  !------------------------------------------------------------------------------
  !> Diffusivity gradient, case 3

  subroutine GetDiffusivityGradient_Case_3(n, x, t, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable diffusivity
    real(RNP), intent(out) :: grad_nu(n,3) !< diffusivity gradient

    real(RNP), allocatable, save :: strain_rate(:,:,:), hessian_v(:,:,:,:)
    real(RNP) :: coeff
    integer   :: i

    !$omp single
    allocate(strain_rate(n,3,3), hessian_v(n,3,3,3))
    !$omp end single

    call GetStrainRate        (n, x, t, strain_rate)
    call GetVelocityHessian   (n, x, t, hessian_v)

    coeff = 4 * nu_1

    !$omp do
    do i = 1, n

      grad_nu(i,1) = coeff * ( strain_rate(i,1,1) * hessian_v(i,1,1,1)   &
                             + strain_rate(i,1,2) * hessian_v(i,1,1,2)   &
                             + strain_rate(i,1,3) * hessian_v(i,1,1,3)   &
                             + strain_rate(i,1,2) * hessian_v(i,1,2,1)   &
                             + strain_rate(i,2,2) * hessian_v(i,1,2,2)   &
                             + strain_rate(i,2,3) * hessian_v(i,1,2,3)   &
                             + strain_rate(i,1,3) * hessian_v(i,1,3,1)   &
                             + strain_rate(i,2,3) * hessian_v(i,1,3,2)   &
                             + strain_rate(i,3,3) * hessian_v(i,1,3,3)   )

      grad_nu(i,2) = coeff * ( strain_rate(i,1,1) * hessian_v(i,1,2,1)   &
                             + strain_rate(i,1,2) * hessian_v(i,1,2,2)   &
                             + strain_rate(i,1,3) * hessian_v(i,1,2,3)   &
                             + strain_rate(i,1,2) * hessian_v(i,2,2,1)   &
                             + strain_rate(i,2,2) * hessian_v(i,2,2,2)   &
                             + strain_rate(i,2,3) * hessian_v(i,2,2,3)   &
                             + strain_rate(i,1,3) * hessian_v(i,2,3,1)   &
                             + strain_rate(i,2,3) * hessian_v(i,2,3,2)   &
                             + strain_rate(i,3,3) * hessian_v(i,2,3,3)   )

      grad_nu(i,3) = coeff * ( strain_rate(i,1,1) * hessian_v(i,1,3,1)   &
                             + strain_rate(i,1,2) * hessian_v(i,1,3,2)   &
                             + strain_rate(i,1,3) * hessian_v(i,1,3,3)   &
                             + strain_rate(i,1,2) * hessian_v(i,2,3,1)   &
                             + strain_rate(i,2,2) * hessian_v(i,2,3,2)   &
                             + strain_rate(i,2,3) * hessian_v(i,2,3,3)   &
                             + strain_rate(i,1,3) * hessian_v(i,3,3,1)   &
                             + strain_rate(i,2,3) * hessian_v(i,3,3,2)   &
                             + strain_rate(i,3,3) * hessian_v(i,3,3,3)   )
    end do

    !$omp barrier
    !$omp master
    deallocate(strain_rate, hessian_v)
    !$omp end master

  end subroutine GetDiffusivityGradient_Case_3

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient, case 4, using exact velocity

  subroutine GetDiffusivityGradient_Case_4(n, x, t, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable diffusivity
    real(RNP), intent(out) :: grad_nu(n,3) !< diffusivity gradient

    real(RNP), allocatable, save :: v(:,:), grad_v(:,:,:)
    integer :: i, l

    !$omp single
    allocate(v(n,3), grad_v(n,3,3))
    !$omp end single

    call GetVelocity(n, x, t, v)
    call GetVelocityGradient(n, x, t, grad_v)

    !$omp do
    do i = 1, n
      do l = 1, 3
        grad_nu(i,l) = nu_1 * 2 * ( grad_v(i,l,1) * v(i,1) &
                                  + grad_v(i,l,2) * v(i,2) &
                                  + grad_v(i,l,3) * v(i,3) )
      end do
    end do

    !$omp barrier
    !$omp master
    deallocate(v, grad_v)
    !$omp end master

  end subroutine GetDiffusivityGradient_Case_4

  !-----------------------------------------------------------------------------
  !> MMS source

  subroutine Get_MMS_Source(problem, n, x, t, f)
    class(FlowProblem_VariableViscosity), intent(in)  :: problem
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: f(n,3)    !< MMS source

    real(RNP), allocatable, save :: dt_u(:,:), nonlin(:,:), grad_p(:,:), &
                                    diff(:,:)
    integer :: i, l, c_stokes

    !$omp single
    allocate(dt_u(n,3), nonlin(n,3), grad_p(n,3), diff(n,3))
    !$omp end single

    call GetVelocityTimeDerivative(n, x, t,          dt_u  )
    call GetNonlinearTerm         (n, x, t,          nonlin)
    call GetPressureGradient      (n, x, t,          grad_p)
    call GetDiffusiveTerm         (problem, n, x, t, diff  )

    if(problem % stokes) then
      c_stokes = 0
    else
      c_stokes = 1
    end if

    do l = 1, 3
      !$omp do
      do i = 1, n
        f(i,l) = dt_u(i,l) + c_stokes * nonlin(i,l) + grad_p(i,l) - diff(i,l)
      end do
    end do

    !$omp barrier
    !$omp master
    deallocate(dt_u, nonlin, grad_p, diff)
    !$omp end master

  end subroutine Get_MMS_Source

  !-----------------------------------------------------------------------------
  !> Pressure gradient

  subroutine GetPressureGradient(n, x, t, grad_p)
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: grad_p(n,3)   !< pressure gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      grad_p(i,1) = k * cos(kxt) * sin(kyt) * sin(kzt)
      grad_p(i,2) = k * sin(kxt) * cos(kyt) * sin(kzt)
      grad_p(i,3) = k * sin(kxt) * sin(kyt) * cos(kzt)
    end do

  end subroutine GetPressureGradient

  !-----------------------------------------------------------------------------
  !> Velocity Gradient
  !>
  !>       grad_v(i,l,m) = ∂v_m/∂x_l
  !>
  !> for mesh point `x(i,:)` and time `t`

  subroutine GetVelocityGradient(n, x, t, grad_v)
    integer,   intent(in)  :: n               !< number of points
    real(RNP), intent(in)  :: x(n,3)          !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: grad_v(n,3,3)   !< velocity gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      grad_v(i,1,1) =  k * cos(kxt) * sin(kzt)                  ! ∂u/∂x
      grad_v(i,2,1) = -k * sin(kyt) * sin(kzt)                  ! ∂u/∂y
      grad_v(i,3,1) =  k * ( sin(kxt) + cos(kyt) ) * cos(kzt)   ! ∂u/∂z

      grad_v(i,1,2) = -k * sin(kxt) * sin(kzt)                  ! ∂v/∂x
      grad_v(i,2,2) =  k * cos(kyt) * sin(kzt)                  ! ∂v/∂y
      grad_v(i,3,2) =  k * ( cos(kxt) + sin(kyt) ) * cos(kzt)   ! ∂v/∂z

      grad_v(i,1,3) = -k * sin(kxt) * cos(kzt)                  ! ∂w/∂x
      grad_v(i,2,3) = -k * sin(kyt) * cos(kzt)                  ! ∂w/∂y
      grad_v(i,3,3) = -k * ( cos(kxt) + cos(kyt) ) * sin(kzt)   ! ∂w/∂z
    end do

  end subroutine GetVelocityGradient

  !-----------------------------------------------------------------------------
  !> Laplacian of velocity components
  !>
  !>       laplacian_v(i,j) = ∇²v_j
  !>
  !> for mesh point `x(i,:)` and time `t`

  subroutine GetVelocityLaplacian(n, x, t, laplacian_v)
    integer,   intent(in)  :: n                !< number of points
    real(RNP), intent(in)  :: x(n,3)           !< mesh points
    real(RNP), intent(in)  :: t                !< time
    real(RNP), intent(out) :: laplacian_v(n,3) !< Laplacian of velocity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)

      laplacian_v(i,1) = -2 * k**2 * (sin(kxt)+cos(kyt)) * sin(kzt)  ! ∇²u(x,t)
      laplacian_v(i,2) = -2 * k**2 * (cos(kxt)+sin(kyt)) * sin(kzt)  ! ∇²v(x,t)
      laplacian_v(i,3) = -2 * k**2 * (cos(kxt)+cos(kyt)) * cos(kzt)  ! ∇²w(x,t)
    end do

  end subroutine GetVelocityLaplacian

  !-----------------------------------------------------------------------------
  !> Hessian of velocity components
  !>
  !>       hessian_v(i,l,m,j) = ∂²v_j /(∂x_l ∂x_m)
  !>
  !> for mesh point `x(i,:)` and time `t`

  subroutine GetVelocityHessian(n, x, t, hessian_v)
    integer,   intent(in)  :: n                  !< number of points
    real(RNP), intent(in)  :: x(n,3)             !< mesh points
    real(RNP), intent(in)  :: t                  !< time
    real(RNP), intent(out) :: hessian_v(n,3,3,3) !< hessian of velocity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)

      ! hessian for velocity component u
      hessian_v(i,1,1,1) = -k**2 * sin(kxt) * sin(kzt)
      hessian_v(i,2,2,1) = -k**2 * cos(kyt) * sin(kzt)
      hessian_v(i,3,3,1) = -k**2 * (sin(kxt) + cos(kyt)) * sin(kzt)
      hessian_v(i,2,3,1) = -k**2 * sin(kyt) * cos(kzt)
      hessian_v(i,1,3,1) =  k**2 * cos(kxt) * cos(kzt)
      hessian_v(i,1,2,1) = 0
      hessian_v(i,3,2,1) = hessian_v(i,2,3,1)
      hessian_v(i,3,1,1) = hessian_v(i,1,3,1)
      hessian_v(i,2,1,1) = hessian_v(i,1,2,1)

      ! hessian for velocity component v
      hessian_v(i,1,1,2) = -k**2 * cos(kxt) * sin(kzt)
      hessian_v(i,2,2,2) = -k**2 * sin(kyt) * sin(kzt)
      hessian_v(i,3,3,2) = -k**2 * (cos(kxt) + sin(kyt)) * sin(kzt)
      hessian_v(i,2,3,2) =  k**2 * cos(kyt) * cos(kzt)
      hessian_v(i,1,3,2) = -k**2 * sin(kxt) * cos(kzt)
      hessian_v(i,1,2,2) = 0
      hessian_v(i,3,2,2) = hessian_v(i,2,3,2)
      hessian_v(i,3,1,2) = hessian_v(i,1,3,2)
      hessian_v(i,2,1,2) = hessian_v(i,1,2,2)

      ! hessian for velocity component w
      hessian_v(i,1,1,3) = -k**2 * cos(kxt) * cos(kzt)
      hessian_v(i,2,2,3) = -k**2 * cos(kyt) * cos(kzt)
      hessian_v(i,3,3,3) = -k**2 * (cos(kxt) + cos(kyt)) * cos(kzt)
      hessian_v(i,2,3,3) =  k**2 * sin(kyt) * sin(kzt)
      hessian_v(i,1,3,3) =  k**2 * sin(kxt) * sin(kzt)
      hessian_v(i,1,2,3) = 0
      hessian_v(i,3,2,3) = hessian_v(i,2,3,3)
      hessian_v(i,3,1,3) = hessian_v(i,1,3,3)
      hessian_v(i,2,1,3) = hessian_v(i,1,2,3)

    end do

  end subroutine GetVelocityHessian

  !-----------------------------------------------------------------------------
  !> Diffusion Term

  subroutine GetDiffusiveTerm(problem, n, x, t, diff)
    class(FlowProblem_VariableViscosity), intent(in) :: problem
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: diff(n,3)     !< diffusive term

    real(RNP), allocatable, save :: v(:,:), laplacian_v(:,:), strain_rate(:,:,:)
    real(RNP), allocatable, save :: nu(:,:), grad_nu(:,:)
    real(RNP) :: tmp
    integer   :: i, l, m

    !$omp single
    allocate(v(n,3), laplacian_v(n,3), strain_rate(n,3,3))
    allocate(nu(n,3), grad_nu(n,3))
    !$omp end single

    call GetVelocity          (n, x, t, v)
    call GetVelocityLaplacian (n, x, t, laplacian_v)
    call GetStrainRate        (n, x, t, strain_rate)

    associate(nu_0 => problem % nu_0,  nu_1 => problem % nu_1)
      select case(problem % test_case)
        case(1)
          call GetDiffusivity_Case_1(n, x, nu_0, nu_1, nu)
          call GetDiffusivityGradient_Case_1(n, x, nu_1, grad_nu)
        case(2)
          call GetDiffusivity_Case_2(n, x, t, nu_0, nu_1, nu)
          call GetDiffusivityGradient_Case_2(n, x, t, nu_1, grad_nu)
        case(3)
          call GetDiffusivity_Case_3(n, x, t, nu_0, nu_1, nu)
          call GetDiffusivityGradient_Case_3(n, x, t, nu_1, grad_nu)
        case(4)
          call GetDiffusivity_Case_4(n, nu_0, nu_1, v, nu)
          call GetDiffusivityGradient_Case_4(n, x, t, nu_1, grad_nu)
      end select
    end associate

    ! diffusive term: ∇(nu 2 S) = ∇nu · 2 S + nu ∇²v
    do l = 1, 3
      !$omp do
      do i = 1, n
        tmp = ZERO
        do m = 1, 3
          tmp = tmp + grad_nu(i,m) * 2 * strain_rate(i,m,l)
        end do
        diff(i,l) = tmp + nu(i,1) * laplacian_v(i,l)
      end do
    end do

    !$omp barrier
    !$omp master
    deallocate(v, laplacian_v, strain_rate, nu, grad_nu)
    !$omp end master

  end subroutine GetDiffusiveTerm

  !-----------------------------------------------------------------------------
  !> Nonlinear Term

  subroutine GetNonlinearTerm(n, x, t, nonlin)
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: nonlin(n,3)   !< nonlinear term

    real(RNP), allocatable, save :: grad_v(:,:,:), v(:,:)
    real(RNP) :: tmp
    integer   :: i, l, m

    !$omp single
    allocate(grad_v(n,3,3), v(n,3))
    !$omp end single

    call GetVelocity        (n, x, t, v      )
    call GetVelocityGradient(n, x, t, grad_v )

    do m = 1, 3
      !$omp do
      do i = 1, n
        tmp = ZERO
        do l = 1, 3
          tmp = tmp + v(i,l) * grad_v(i,l,m)
        end do
        nonlin(i,m) = tmp
      end do
    end do

    !$omp barrier
    !$omp master
    deallocate(grad_v, v)
    !$omp end master

  end subroutine GetNonlinearTerm

  !-----------------------------------------------------------------------------
  !> Strain rate tensor

  subroutine GetStrainRate(n, x, t, strain_rate)
    integer,   intent(in)  :: n                    !< number of points
    real(RNP), intent(in)  :: x(n,3)               !< mesh points
    real(RNP), intent(in)  :: t                    !< time
    real(RNP), intent(out) :: strain_rate(n,3,3)   !< strain rate tensor

    real(RNP), allocatable, save :: grad_v(:,:,:)
    integer   :: i, l, m

    !$omp single
    allocate(grad_v(n,3,3))
    !$omp end single

    call GetVelocityGradient(n, x, t, grad_v)

    do m = 1, 3
    do l = 1, 3
      !$omp do
      do i = 1, n
        ! symmetric part of velocity gradient tensor: S = 1/2 * (∇v + (∇v)ᵀ)
        strain_rate(i,l,m) = HALF * (grad_v(i,l,m) + grad_v(i,m,l))
      end do
    end do
    end do

    !$omp barrier
    !$omp master
    deallocate(grad_v)
    !$omp end master

  end subroutine GetStrainRate

!===============================================================================

end module ISP_Flow_Problem__Variable_Viscosity

