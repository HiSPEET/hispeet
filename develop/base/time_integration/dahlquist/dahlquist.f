!> summary:  Stability and accuracy analysis based on the Dahlquist equation
!> author:   Joerg Stiller and Martina Grotteschi
!> date:     2020/07/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Stability and accuracy analysis based on the Dahlquist equation
!>
!> The program solves the Dahlquist test equation
!>
!>     ∂z/∂t = λ z,      with λ ∈ c, Re(λ) ∈ R⁻
!>
!> The domain of stability is defined as the set S = { z ∈ c : |z| ≤ 1 }.
!>
!> The following methods are available for stability analysis:
!>
!>   -  Euler (Eu)
!>   -  Trapezoidal rule (TR)
!>   -  Runge-Kutta (RK)
!>   -  SDC with Eu/TR/RK predictor and Eu/RK corrector using either equidistant
!>      or Lobatto collocation points
!>
!> All methods are available as implicit, explicit and semi-implicit (IMEX)
!> variants. With SDC the degree of implicitness can be chosen independently
!> for the predictor and for the corrector.
!>
!> In order to solve the test equation, a single time step of length 1 is
!> executed.
!>
!>#### Usage
!>
!> The program is started from command line with the case identifier as one
!> optional argument, e.g.
!>
!>     ./dahlquist
!>
!> The input parameters are specified in namelist records contained in a single
!> file named according to the case identifier with suffix `.prm`. In the above
!> example the program would search for an input file named `dahlquist.prm`.
!>
!>   *  Re(λ) outputs,
!>   *  Im(λ) outputs,
!>   *  a,
!>   *  e, where `e` is the error.
!>
!> The results can be evaluated using the Python scripts `stability.py` and
!> `accuracy.py`
!===============================================================================

program Dahlquist
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use DQ__Time_Integrator
  use DQ__Time_Integrator__Euler
  use DQ__Time_Integrator__RK
  use DQ__Time_Integrator__TR
  use DQ__SDc__Method
  use DQ__SDc__Method__Euler
  use DQ__SDc__Method__RK

  implicit none

  integer :: time_method = 1 ! standalone or predictor method
                             ! 1  Euler
                             ! 2  Trapezoidal rule
                             ! 3  Runge-Kutta

  integer :: sdc_method  = 0 ! SDC method
                             ! 0  none (use standalone time integrator)
                             ! 1  based on Euler
                             ! 3  based on Runge-Kutta (default)

  namelist /input/ time_method, sdc_method

  ! standalone time-integrator or predictor
  class(DQ_TimeIntegrator), allocatable :: tint
  type(DQ_TimeIntegrator_Options_Euler) :: tint_opt_eu ! options for Euler
  type(DQ_TimeIntegrator_Options_TR)    :: tint_opt_tr ! options for TR
  type(DQ_TimeIntegrator_Options_RK)    :: tint_opt_rk ! options for RK

  namelist /input/ tint_opt_eu, tint_opt_tr, tint_opt_rk

  ! SDC
  class(DQ_SDc_Method), allocatable :: sdc
  type(DQ_SDc_Options_Euler)        :: sdc_opt_eu ! options for Euler-based SDC
  type(DQ_SDc_Options_RK)           :: sdc_opt_rk ! options for RK-based SDC

  namelist /input/ sdc_opt_eu, sdc_opt_rk

  real(RNP) :: c_min =   0  ! min CFL number
  real(RNP) :: c_max =  10  ! max CFL number
  real(RNP) :: d_min =  -5  ! min diffusion number
  real(RNP) :: d_max =  10  ! max diffusion number
  integer   :: nc    =  80  ! num intervals along imaginary axis
  integer   :: nd    =  80  ! num intervals along real axis

  namelist /input/ c_min, c_max, d_min, d_max, nc, nd

  real(RNP)   , allocatable :: c(:) ! convective CFL numbers (0:nc)
  real(RNP)   , allocatable :: d(:) ! diffusive  CFL numbers (0:nd)

  complex(RNP), allocatable :: lambda(:,:)  ! λ = -d + ic = λᵣ + iλᵢ
  real(RNP)   , allocatable :: a(:,:)       ! amplification: a = |z|
  real(RNP)   , allocatable :: e(:,:)       ! error        : e = |z -z_ex|

  character(len=80) :: input_file = 'dahlquist'
  complex(RNP) :: z
  real(RNP)    :: delta_c, delta_d
  integer      :: i, j, io

  ! read options and parameters ................................................

  input_file = trim(input_file) // '.prm'

  open(newunit=io, file=input_file)
  read(io, nml=input)
  close(io)

  ! initialize parameters ......................................................

  allocate(c(0:nc), d(0:nd), lambda(0:nd,0:nc), a(0:nd,0:nc), e(0:nd,0:nc))

  ! cfl numbers
  delta_c = (c_max - c_min) / nc
  do j = 0, nc
    c(j) = c_min + j * delta_c
  end do

  ! diffusion numbers
  delta_d = (d_max - d_min) / nd
  do i = 0, nd
    d(i) = d_min + i * delta_d
  end do

  ! λ = -d + ic
  do j = 0, nc
  do i = 0, nd
    lambda(i,j) =  -d(nd - i) + (ZERO,ONE) * c(j)
  end do
  end do

  ! initialize time-integration method .........................................

  select case(sdc_method)

  case(0) ! standalone time integrator
    select case(time_method)
    case(1)
      tint = DQ_TimeIntegrator_Euler(tint_opt_eu)
    case(2)
      tint = DQ_TimeIntegrator_TR(tint_opt_tr)
    case(3)
      tint = DQ_TimeIntegrator_RK(tint_opt_rk)
    end select

  case(1) ! SDC based on Euler
    select case(time_method)
    case(1)
      sdc = DQ_SDC_Method_Euler(tint_opt_eu, sdc_opt_eu)
    case(2)
      sdc = DQ_SDC_Method_Euler(tint_opt_tr, sdc_opt_eu)
    case(3)
      sdc = DQ_SDC_Method_Euler(tint_opt_rk, sdc_opt_eu)
    end select

  case default ! SDC based on Runge-Kutta
    select case(time_method)
    case(1)
      sdc = DQ_SDC_Method_RK(tint_opt_eu, sdc_opt_rk)
    case(2)
      sdc = DQ_SDC_Method_RK(tint_opt_tr, sdc_opt_rk)
    case(3)
      sdc = DQ_SDC_Method_RK(tint_opt_rk, sdc_opt_rk)
    end select

  end select

  ! show settings
  select case(sdc_method)
  case(0)
    call tint % Show()
  case default
    call sdc % Show()
  end select

  ! integration with dt = 1 ....................................................

  do j = 0, nc
  do i = 0, nd

    z = (1,0)

    select case(sdc_method)
    case(0) ! standalone time-integrator
      call tint % TimeStep(lambda(i,j), ONE, z)
    case default ! SDC method
      call sdc % TimeStep(lambda(i,j), ONE, z)
    end select

    a(i,j) = abs(z)                     ! amplification
    e(i,j) = abs(z - exp(lambda(i,j)))  ! error

  end do
  end do

  ! save results ...............................................................

  ! x = Re(λ) = -d
  open(newunit=io, file = 'lambda_re.dat')
  do i = 0, nd
    write(io, '(ES12.5,1X)', advance='NO') lambda(i,0) % re
  end do
  close(io)

  ! y = Im(λ) = c
  open(newunit=io, file = 'lambda_im.dat')
  do j = 0, nc
    write(io, '(ES12.5,1X)', advance='NO') lambda(0,j) % im
  end do
  close(io)

  ! amplification factors
  open(newunit=io, file = 'amplification.dat')
  do j = 0, nc
    do i = 0, nd
      write(io, '(ES12.5,1X)', advance='NO') a(i,j)
    end do
    write(io,*)
  end do
  close(io)

  ! error
  open(newunit=io, file = 'error.dat')
  do j = 0, nc
    do i = 0, nd
    ! write(io, '(ES17.10,1X)', advance='NO') e(i,j)
      write(io, '(ES12.5,1X)', advance='NO') e(i,j)
    end do
    write(io,*)
  end do
  close(io)

  !=============================================================================

end program Dahlquist
