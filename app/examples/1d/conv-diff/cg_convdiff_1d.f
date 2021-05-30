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
!>   *  SDC based on IMEX Euler
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
  use CG_Element_Operators__1D
  use CG_Utilities_1D
  use Spectral_Deferred_Correction
  use Harmonic_Wave_Package

  use CG_ConvDiff_1D__IMEX_Euler
  use CG_ConvDiff_1D__IMEX_BDF2
  use CG_ConvDiff_1D__IMEX_BDF3
  use CG_ConvDiff_1D__IMEX_RK
  use CG_ConvDiff_1D__IMEX_Euler_SDC
  use CG_ConvDiff_1D__RK_SDC
  use CG_ConvDiff_1D__IMEX_TR_SDC

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! problem ....................................................................

  ! case identifier: first argument or, if absent, name of inving command
  character(len=200) :: conv_diff_case
  character(len=204) :: input_file   ! = trim(conv_diff_case) // '.prm'
  character(len=204) :: result_file  ! = trim(conv_diff_case) // '.dat'


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
  integer   :: po = 16         ! polynomial order
  integer   :: ne = 10         ! number of elements
  real(RNP) :: dx              ! element width

  namelist /discretization_parameters/ po, ne

  ! time                                               priority:
  real(RNP) :: cfl    = -1     ! CFL number              1  if > 0 and v ≠ 0
  real(RNP) :: dt     = -1     ! time step size          2  if > 0
  integer   :: nt     = -1     ! number of time steps    3  upper limit if > 0
  integer   :: method =  1     ! time stepping method
  ! available methods
  ! 1   IMEX Euler
  ! 2   IMEX BDF2
  ! 3   IMEX BDF3
  ! 4   IMEX Runge-Kutta
  ! 5   SDC based on IMEX Euler
  logical :: flying_start = .false. ! use exact solution for t < 0

  namelist /discretization_parameters/ cfl, dt, nt, method, flying_start


  ! operators ..................................................................

  type(CG_ElementOperators_1D) :: eop          ! element operators

  type(IMEX_RK_Method1D)      :: imex_rk      ! IMEX Runge-Kutta method

  type(SDC_Method)            :: sdc          ! Euler SDC method
  type(SDC_Options)           :: sdc_opt      ! Euler SDC options
  namelist /sdc_parameters/      sdc_opt

  type(RK_SDC_Method1D)       :: rk_sdc       ! IMEX RK method

  real(RNP), allocatable      :: M(:,:)       ! global mass matrix

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

  logical   :: exists, periodic
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

    select case (method)
    case(4)
      call Init_IMEX_RK(imex_rk, po, ne, io)
    case(5)
      call Init_Euler_SDC(sdc, io)
    case(6)
      call Init_RK_SDC(rk_sdc, po, ne, io)
    case(7)
      call Init_TR_SDC(sdc, io)
    end select
    rewind(io)
    call wave % New(input_file)
    close(io)
  else
    call Error('ConvDiff_CG_SEM_1D', &
               'case file "' // trim(input_file) // '" not found')
  end if
  write(*,'(/,A,1X,I0)') 'method =', method

  periodic = all(bc == 'P')

  ! workspace ..................................................................

  allocate(M    (0:po,ne))
  allocate(x    (0:po,ne)) ! since mold is dangerous with lower bound of 0
  allocate(w    (0:po,ne)) !
  allocate(u    (0:po,ne)) !
  allocate(u0   (0:po,ne)) !
  allocate(u1   (0:po,ne)) !
  allocate(u2   (0:po,ne)) !
  allocate(u_ex (0:po,ne)) !
  allocate(err  (0:po,ne)) !

  ! mesh and operators .........................................................

  eop = CG_ElementOperators_1D(po)
  call GetMeshPoints(eop, ZERO, ONE, dx, x)
  call GetMassMatrix(eop, dx, M, periodic)
  call GetPointWeights(w, periodic)

  ! initial conditions .........................................................

  t = 0
  call wave % GetAmplitude(v, nu, x, t       , u0)
  call wave % GetAmplitude(v, nu, x, t -   dt, u1)
  call wave % GetAmplitude(v, nu, x, t - 2*dt, u2)

  ! time integration ...........................................................

  call SetTimeIntegrationParameters(po, dx, v, t_end, cfl, dt, nt)
  write(*,'(/,A,1X,I0,/)') 'nt =', nt

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

      case(5)
        call IMEX_Euler_SDC(sdc, eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)

      case(6)
        call rk_sdc % TimeStep(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)

      case(7)
        call IMEX_TR_SDC(sdc, eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)


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

  write(*,'(9(A8,7X))') '#      t', 'dt', 'cfl', 'err_max', 'err_2'
  write(*,'(5(ES12.5,3X),A)') t, dt, cfl, err_max, err_2, '#last#'

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

  subroutine Init_IMEX_RK(imex_rk, po, ne, io)
    type(IMEX_RK_Method1D), intent(inout) :: imex_rk !< IMEX RK method
    integer, intent(in) :: po !< polynomial order
    integer, intent(in) :: ne !< number of elements
    integer, intent(in) :: io !< unit number of input file

    integer :: ns     = 3 ! number of stages
    integer :: method = 1 ! RK method, if several with `ns`stages exist
    logical :: show   = .false. ! print IMEX RK properties and coefficients

    namelist /imex_rk_parameters/ ns, method, show

    read(io, nml=imex_rk_parameters)
    imex_rk = IMEX_RK_Method1D(po, ne, ns, method)

    if (show) then
      call imex_rk % Show()
    end if

  end subroutine Init_IMEX_RK

  !-----------------------------------------------------------------------------
  !> Initialization of the Euler-SDC method

  subroutine Init_Euler_SDC(sdc, io)
    type(SDC_Method), intent(inout) :: sdc !< SDC method
    integer, intent(in) :: io !< unit number of input file

    integer   :: n_sub     =  1      ! number of subintervals (M)
    integer   :: n_sweep   =  0      ! max num correction sweeps (K)
    integer   :: point_set =  1      ! equidistant (1) or GLL (2) points
    namelist /euler_sdc_parameters/ n_sub, n_sweep, point_set

    type(SDC_Options) :: opt ! Euler-SDC options

    read(io, nml=euler_sdc_parameters)

    opt = SDC_Options( n_sub     = n_sub,       &
                       n_sweep   = n_sweep,     &
                       point_set = point_set    )

    sdc = SDC_Method(opt)

  end subroutine Init_Euler_SDC


  !-----------------------------------------------------------------------------
  !> Initialization of the IMEX-RK-SDC method

  subroutine Init_RK_SDC(rk_sdc, po, ne, io)
    type(RK_SDC_Method1D), intent(inout) :: rk_sdc !< RK-SDC method
    integer, intent(in) :: po !< polynomial order
    integer, intent(in) :: ne !< number of elements
    integer, intent(in) :: io !< unit number of input file

    ! SDC parameters
    integer   :: n_sub      =  1    ! number of subintervals (M)
    integer   :: n_sweep    =  0    ! max num correction sweeps (K)
    integer   :: point_set  =  2    ! equidistant (1) or GLL (2) points
    namelist /rk_sdc_parameters/ n_sub, n_sweep, point_set

    ! IMEX Runge-Kutta parameters
    integer :: n_stage = 3       ! number of stages
    integer :: method  = 1       ! RK method, if several with `ns`stages exist
    logical :: show    = .false. ! print IMEX RK properties and coefficients
    namelist /rk_sdc_parameters/ n_stage, method, show

    type(SDC_Options) :: sdc_opt ! RK-SDC options

    read(io, nml=rk_sdc_parameters)

    sdc_opt = SDC_Options( n_sub     = n_sub,     &
                           n_sweep   = n_sweep,   &
                           point_set = point_set  )

    rk_sdc = RK_SDC_Method1D(po, ne, sdc_opt, n_stage, method)

    if (show) then
      call rk_sdc % imex_rk % Show()
    end if

  end subroutine Init_RK_SDC

 !-----------------------------------------------------------------------------
 !> Initialization of the IMEX-TR-SDC method


  subroutine Init_TR_SDC(sdc, io)
    type(SDC_Method), intent(inout) :: sdc !< SDC method
    integer, intent(in) :: io !< unit number of input file

    ! TR-SDC parameters
    integer   :: n_sub      =  1      ! number of subintervals (M)
    integer   :: n_sweep    =  0      ! max num correction sweeps (K)
    integer   :: point_set  =  2      ! equidistant (1) or GLL (2) points

    logical :: show   = .false. ! print IMEX TR properties and coefficients

    namelist /tr_sdc_parameters/ n_sub, n_sweep, point_set, show

    type(SDC_Options) :: opt ! SDC options

    read(io, nml=tr_sdc_parameters)

    opt = SDC_Options( n_sub     = n_sub,     &
                       n_sweep   = n_sweep,   &
                       point_set = point_set  )

    sdc = SDC_Method(opt)

    !if (show) then
     ! call sdc % Write()
    !end if


  end subroutine Init_TR_SDC


!==============================================================================

end program CG_ConvDiff_1D
