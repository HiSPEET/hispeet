!> summary:  Spectral element solver for linear convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral element solver for linear convection-diffusion equation
!>
!> The program solves the convection-diffusion equation
!>
!>     ∂u/∂t + v ∂u/∂x = nu ∂²u/∂x²
!>
!> in the domain Ω = (0,1) for 0 < t ≤ t_end. Either periodicity `bc='P'`, or a
!> combination of Dirichlet (`'D'`) and Neumann (`'N'`) boundary conditions can
!> be imposed.
!>
!> The initial and boundary values are derived from the exact solution, which is
!> given in the form of a harmonic wave package
!>
!>     u(x,t) = ∑ᵢ aᵢ sin[2πkᵢ/lw (x - sᵢ - vt)] exp[-(2πkᵢ)² nu t]
!>
!> Spatial discretization is based on the continuous Galerkin spectral-element
!> using nodal base functions and GLL quadrature.
!>
!> The following methods are available for time integration:
!>
!>   *  IMEX Euler
!>   *  IMEX BDF2
!>   *  IMEX BDF3
!>   *  IMEX Runge-Kutta
!>
!>#### Usage
!>
!> The program is started from command line with the case identifier as one
!> optional argument, e.g.
!>
!>     cg_convdiff_1d my_test
!>
!> In this example, the case identifier is `my_test`. If not specified, the
!> identifier is set to the name of the execuatble, i.e. `cg_convdiff_1d`.
!>
!> The input parameters are specified in namelist records contained in a single
!> file named according to the case identifier with suffix `.prm`. In the above
!> example the program would search for an input file named `my_test.prm`.
!>
!> The input file must provide the following namelist records
!>
!>   *  problem_parameters
!>   *  wave_package_dimensions
!>   *  wave_parameters
!>   *  discretization_parameters
!>
!> and, additionally, if an IMEX Runge-Kutta method is used
!>
!>   *  imex_runge_kutta
!>
!> For the description of the input parameters see below and corresponding
!> sections of modules Harmonic_Wave_Package and IMEX_Runge_Kutta_Method.
!> An example is given in the default input file `cg_convdiff_1d.prm`.
!>
!> Upon successful execution the program writes an output file named according
!> to the case identifier with suffix `.dat`. It contains
!>
!>   *  the mesh points `x`,
!>   *  the numerical solution `u`,
!>   *  the exact solution `u_ex`, and
!>   *  the error `err`
!>
!> after the final time step.
!>
!===============================================================================

program CG_ConvDiff_1D
  use Kind_Parameters,   only: RNP, IXL
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error
  use CG_Element_Operators_1D
  use CG_Utilities_1D
  use Harmonic_Wave_Package

  use CG_ConvDiff_1D__IMEX_Euler
  use CG_ConvDiff_1D__IMEX_BDF2
  use CG_ConvDiff_1D__IMEX_BDF3
  use CG_ConvDiff_1D__IMEX_RK

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! problem ....................................................................

  ! case identifier: first argument or, if absent, name of inving command
  character(len=80) :: conv_diff_case
  character(len=84) :: input_file   ! = trim(conv_diff_case) // '.prm'
  character(len=84) :: result_file  ! = trim(conv_diff_case) // '.dat'


  ! problem parameters
  real(RNP) :: v     = 1         ! convection velocity
  real(RNP) :: nu    = 0.001     ! diffusivity
  real(RNP) :: t_end = 1         ! integration time
  character :: bc(2) = ['D','N'] ! left/right BC ('D': Dirichlet, 'N': Neumann)

  namelist /problem_parameters/ v, nu, t_end, bc

  ! wave package representing the exact solution, initialized through namelists
  ! `wave_package_dimensions` and `wave_parameters`, provided in the case file
  type(HarmonicWavePackage) :: wave

  ! discretization .............................................................

  ! space
  integer   :: po = 16       ! polynomial order
  integer   :: ne = 10       ! number of elements
  real(RNP) :: dx            ! element width

  namelist /discretization_parameters/ po, ne

  ! time                                               priority:
  real(RNP) :: cfl    = -1   ! CFL number              1  if > 0 and v ≠ 0
  real(RNP) :: dt     = -1   ! time step size          2  if > 0
  integer   :: nt     = -1   ! number of time steps    3  upper limit if > 0
  integer   :: method =  1   ! time stepping method
  ! available methods
  ! 1   IMEX Euler
  ! 2   IMEX BDF2
  ! 3   IMEX BDF3
  ! 4   IMEX Runge-Kutta
  logical   :: flying_start = .false. ! use exact solution for t < 0

  namelist /discretization_parameters/ cfl, dt, nt, method, flying_start

  ! operators ..................................................................

  type(CG_ElementOperators1D) :: eop     ! element operators
  real(RNP), allocatable      :: M(:,:)  ! global mass matrix

  type(ConvDiff_IMEX_RK)      :: imex_rk ! IMEX Runge-Kutta method


  ! variables ..................................................................

  real(RNP)              :: t           ! time
  real(RNP), allocatable :: x    (:,:)  ! mesh points
  real(RNP), allocatable :: w    (:,:)  ! point weights
  real(RNP), allocatable :: u    (:,:)  ! solution at t = t₀ +  ∆t
  real(RNP), allocatable :: u0   (:,:)  ! solution at t = t₀
  real(RNP), allocatable :: u1   (:,:)  ! solution at t = t₀ -  ∆t
  real(RNP), allocatable :: u2   (:,:)  ! solution at t = t₀ - 2∆t
  real(RNP), allocatable :: u_ex (:,:)  ! exact solution
  real(RNP), allocatable :: err  (:,:)  ! error

  ! auxiliary ..................................................................

  logical   :: exists
  integer   :: io, stat
  integer   :: i, e
  real(RNP) :: err_max, err_2

  !-----------------------------------------------------------------------------
  ! Intitialization

  ! greeting
  write(*,'(/,A)') 'CG-SEM for the 1D convection-diffusion equation'

  ! identify case ..............................................................

  call get_command_argument(1, conv_diff_case, status=stat)
  if (stat /= 0 .or. len_trim(conv_diff_case) == 0) then
    call get_command_argument(0, conv_diff_case, status=stat)
  end if
  input_file  = trim(conv_diff_case) // '.prm'
  result_file = trim(conv_diff_case) // '.dat'

  ! load parameters ............................................................

  inquire(file=input_file, exist=exists)
  if (exists) then
    open(newunit=io, file=input_file)
    read(io, nml=problem_parameters)
    read(io, nml=discretization_parameters)
    rewind(io)
    call wave % New(input_file)
    close(io)
  else
    call Error('ConvDiff_CG_SEM_1D', &
               'case file "' // trim(input_file) // '" not found')
  end if

  ! workspace ..................................................................

  allocate(M(0:po,ne))
  allocate(x,    mold = M) ! join
  allocate(w,    mold = M) ! these
  allocate(u,    mold = M) ! once
  allocate(u0,   mold = M) ! supported
  allocate(u1,   mold = M) ! by
  allocate(u2,   mold = M) ! ..
  allocate(u_ex, mold = M) ! ..
  allocate(err,  mold = M) ! pgfortran

  ! mesh and operators .........................................................

  call eop % New(po)
  call eop % BuildInteriorEigensystem()
  call GetMeshPoints(eop, ZERO, ONE, dx, x)
  call GetMassMatrix(eop, dx, bc, M)
  call GetPointWeights(bc, w)

  ! initial conditions .........................................................

  t = 0
  call wave % GetAmplitude(v, nu, x, t       , u0)
  call wave % GetAmplitude(v, nu, x, t -   dt, u1)
  call wave % GetAmplitude(v, nu, x, t - 2*dt, u2)

  ! time integration ...........................................................

  call SetTimeIntegrationParameters(po, dx, v, t_end, cfl, dt, nt)
  write(*,'(/,A,1X,I0,/)')    'nt =', nt

  select case(method)
  case(4)
    call Init_IMEX_RK(imex_rk, input_file, po, ne)
  end select

  !-----------------------------------------------------------------------------
  ! Time integration

  if (nt > 0) then

    do i = 1, nt

      select case(method)

      case(1)
        call IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)

      case(2)
        if (i == 1 .and. .not. flying_start) then
          call IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)
        else
          call IMEX_BDF2(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u1, u)
        end if

      case(3)
        if (i == 1 .and. .not. flying_start) then
          call IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)
        else if (i == 2 .and. .not. flying_start) then
          call IMEX_BDF2(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u1, u)
        else
          call IMEX_BDF3(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u1, u2, u)
        end if

      case(4)
        call imex_rk % TimeStep(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)

      end select

      t  = t + dt
      u2 = u1
      u1 = u0
      u0 = u

    end do

  else
    u = u0
  end if

  !-----------------------------------------------------------------------------
  ! Evaluation

  ! exact solution
  call wave % GetAmplitude(v, nu, x, t, u_ex)

  ! error
  err = u - u_ex
  err_max = maxval(abs(err))
  err_2   = sqrt(sum(w * M * err**2))

  write(*,'(9(A8,7X))') 'dt', 'cfl', 'err_max', 'err_2'
  write(*,'(9(ES12.5,3X))') dt, cfl, err_max, err_2

  ! save results
  open(newunit=io, file=result_file)
  write(io,'(A1,A11,9(1X,A12))') '#','x', 'u', 'u_ex', 'err'
  do e = 1, ne
  do i = 0, po
    write(io,'(10(ES12.5,1X))') x(i,e), u(i,e), u_ex(i,e), err(i,e)
  end do
  end do
  close(io)

contains

  !-----------------------------------------------------------------------------
  !> Set CFL number, time step width and number of time steps

  subroutine SetTimeIntegrationParameters(po, dx, v, t_end, cfl, dt, nt)
    integer,   intent(in)    :: po    !< polynomial order
    real(RNP), intent(in)    :: dx    !< element length
    real(RNP), intent(in)    :: v     !< velocity
    real(RNP), intent(in)    :: t_end !< integration time
    real(RNP), intent(inout) :: cfl   !< CFL number      (admissible → chosen)
    real(RNP), intent(inout) :: dt    !< time step width (admissible → chosen)
    integer,   intent(inout) :: nt    !< number of steps (admissible → chosen)

    real(RNP) :: h

    h = dx / po

    if (cfl > 0 .and. v /= 0) then
      dt = cfl * h / abs(v)
    else if (dt <= 0 .and. nt > 0) then
      dt = t_end / nt
    else
      dt = max(dt, ZERO)
    end if

    if (dt > 0) then
      cfl = dt * abs(v) / h
      if (nt > 0) then
        nt = min(nt, nint(t_end / dt))
      else
        nt = nint(t_end / dt)
      end if
    else
      cfl = 0
      dt  = 0
      nt  = 0
    end if

  end subroutine SetTimeIntegrationParameters

  !-----------------------------------------------------------------------------
  !> Initialization of the IMEX Runge-Kutta method
  !>
  !> The number of stages `s` and, optionally, the scheme identifier `m` are
  !> read from namelist `imex_runge_kutta` provided in `imex_rk_file`.

  subroutine Init_IMEX_RK(imex_rk, imex_rk_file, po, ne)
    class(ConvDiff_IMEX_RK), intent(inout) :: imex_rk !< IMEX RK method
    character(len=*), intent(in) :: imex_rk_file !< input file
    integer, intent(in) :: po !< polynomial order
    integer, intent(in) :: ne !< number of elements

    logical :: opened, exists
    integer :: io
    integer :: s = 3
    integer :: m = 1

    namelist /imex_runge_kutta/ s, m

    inquire(file=imex_rk_file, opened=opened, exist=exists, number=io)

    if (.not. opened) then
      if (exists) then
        open(newunit=io, file=imex_rk_file, action='READ')
      else
        call Error('Init_IMEX_RK',                                            &
                   'input file "' // trim(imex_rk_file) // '" does not exist' )
      end if
    end if

    read(io, nml=imex_runge_kutta)
    call imex_rk % New(s, m, po, ne)

    ! close IO unit if file was closed on entry
    if (.not. opened) then
      close(io)
      end if

  end subroutine Init_IMEX_RK

!==============================================================================

end program CG_ConvDiff_1D
