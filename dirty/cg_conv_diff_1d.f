!> summary:  Spectral element solver for linear convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral element solver for linear convection-diffusion equation
!>
!>#### Problem
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
!>#### Spatial discretization
!>
!> TBD
!>
!>#### Time discretization
!>
!> TBD
!>
!===============================================================================

program CG_ConvDiff_1D
  use Kind_Parameters,   only: RNP, IXL
  use Constants,         only: ONE, ZERO
  use Execution_Control, only: Error
  use CG_Element_Operators_1D
  use CG_Spectral_Element_Utils_1D
  use Harmonic_Wave_Package

  use CG_ConvDiff_1D__IMEX_Euler
  use CG_ConvDiff_1D__IMEX_BDF2
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

  namelist /discretization_parameters/ cfl, dt, nt, method

  ! operators ..................................................................

  type(CG_ElementOperators1D) :: eop    ! element operators
  real(RNP), allocatable      :: M(:,:) ! global mass matrix

  ! variables ..................................................................

  real(RNP)              :: t = 0       ! time
  real(RNP), allocatable :: x    (:,:)  ! mesh points
  real(RNP), allocatable :: w    (:,:)  ! point weights
  real(RNP), allocatable :: u    (:,:)  ! new solution
  real(RNP), allocatable :: u0   (:,:)  ! old solution
  real(RNP), allocatable :: u1   (:,:)  ! solution two steps before
  real(RNP), allocatable :: u_ex (:,:)  ! exact solution
  real(RNP), allocatable :: err  (:,:)  ! error

  ! auxiliary ..................................................................

  logical :: exists
  integer :: io, stat
  integer :: i, e

  !-----------------------------------------------------------------------------
  ! Intitialization

  ! greeting
  write(*,'(/,A)') 'CG-SEM for the 1D convection-diffusion equation'

  ! identify case and input/result files
  call get_command_argument(1, conv_diff_case, status=stat)
  if (stat /= 0 .or. len_trim(conv_diff_case) == 0) then
    call get_command_argument(0, conv_diff_case, status=stat)
  end if
  input_file  = trim(conv_diff_case) // '.prm'
  result_file = trim(conv_diff_case) // '.dat'

  ! load parameters
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

  ! workspace
  allocate(M(0:po,ne))
  allocate(x,    mold = M) ! join
  allocate(w,    mold = M) ! these
  allocate(u,    mold = M) ! once
  allocate(u0,   mold = M) ! supported
  allocate(u1,   mold = M) ! by
  allocate(u_ex, mold = M) ! pgfortran
  allocate(err,  mold = M) !

  ! mesh and operators
  dx = ONE / ne
  call eop % New(po)
  call eop % InitSubstructuring()
  call GetMeshPoints(eop, dx, x)
  call GetMassMatrix(eop, dx, bc, M)
  call GetPointWeights(bc, w)

  ! initial conditions
  call wave % GetAmplitude(v, nu, x, t, u0)

  ! time integration parameters
  if (cfl > 0 .and. v /= 0) then
    dt = cfl * dx / (po * abs(v))
  else if (dt <= 0 .and. nt > 0) then
    dt = t_end / nt
  else
    dt = max(dt, ZERO)
  end if
  if (dt > 0) then
    cfl = dt * abs(v) * po / dx
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

  write(*,'(/,A,ES12.5)') 'dt  =', dt
  if (v /= 0) then
    write(*,'(A,ES12.5)') 'cfl =', cfl
  end if
  write(*,'(A,1X,I0)')    'nt  =', nt

  !-----------------------------------------------------------------------------
  ! Time integration

  if (nt > 0) then

    do i = 1, nt
      select case(method)
      case(1)
        call IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)
      case(2)
        if (i == 1) then
          call IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u)
        else
          call IMEX_BDF2(eop, dx, dt, M, wave, v, nu, bc, x, t, u0, u1, u)
        end if
      end select
      t  = t + dt
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
  write(*,'(/,A)') 'error'
  write(*,'(2X,A,ES12.5)') 'err_max =', maxval(abs(err))
  write(*,'(2X,A,ES12.5)') 'err_2   =', sqrt(sum(w * M * err**2))

  ! save results
  open(newunit=io, file=result_file)
  write(io,'(10(A12,1X))') '# x', 'u', 'u_ex', 'err'
  do e = 1, ne
  do i = 0, po
    write(io,'(10(ES12.5,1X))') x(i,e), u(i,e), u_ex(i,e), err(i,e)
  end do
  end do
  close(io)


!==============================================================================

end program CG_ConvDiff_1D
