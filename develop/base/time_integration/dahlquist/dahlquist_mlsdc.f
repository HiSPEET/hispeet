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

!> summary:  Dahlquist Problem MLSDC
!> author:   Erik Pfister
!> date:     2024/04/24
!===============================================================================

program Dahlquist_MLSDC
  use Kind_Parameters
  use Constants
  use Execution_Control

  use DQ__Time_Integrator
  use DQ__Time_Integrator__Euler
  use DQ__Time_Integrator__ISD
  use DQ__SDC__Method
  use DQ__MLSDC__Corrector__Euler
  use DQ__MLSDC__Corrector__Euler_BEPar
  use DQ__MLSDC__Corrector__ISD
  use DQ__MLSDC__1D
  use DQ__MLSDC__Level__1D
  use DQ__MLSDC__Variable__1D
  use DQ__MLSDC__Cascade__1D
  use DQ__MLSDC__FMG__1D
  use DQ__MLSDC__V_Cycle__1D

  implicit none

  ! declarations: control ......................................................

  character(len=80) :: case_name = 'dahlquist_mlsdc'
  character(len=80) :: case_file

  ! declarations: discretization ...............................................

  integer, parameter :: max_n_level = 20 ! upper bound for number of levels

  real(RNP) :: dt_slab  = 1.0   ! thickness of one time slab
  integer   :: n_level  = 2     ! number of space-time levels
  integer   :: n_cycle  = 2    ! number of v-cycles
  integer   :: n_coarse = 2    ! number of coarse sweeps

  namelist/discretization_prm/ dt_slab, n_level, n_cycle, n_coarse

  integer :: p_time  (max_n_level)  = -1 ! polynomial degree of time step
  integer :: n_time  (max_n_level)  = -1 ! number of time steps in one slice
  integer :: n_stage_p(max_n_level) = -1
  integer :: n_stage_c(max_n_level) = -1

  namelist/discretization_prm/ p_time, n_time, n_stage_p, n_stage_c

  ! time integration -- still needs to be configured
  class(DQ_TimeIntegrator_Options), allocatable :: opt_pre
  class(DQ_SDC_Options)           , allocatable :: opt_sdc

  logical :: incremental = .true.
  integer :: mg_start    = 1
  integer :: time_method = 3
  integer :: sdc_method  = 2

  namelist/time_integration_prm/ incremental, mg_start, time_method, sdc_method

  type(DQ_TimeIntegrator_Options_Euler)        :: opt_pre_euler
  type(DQ_TimeIntegrator_Options_ISD)          :: opt_pre_isd
  type(DQ_MLSDC_Corrector_Options_Euler)       :: opt_sdc_euler
  type(DQ_MLSDC_Corrector_Options_Euler_BEPar) :: opt_sdc_euler_bepar
  type(DQ_MLSDC_Corrector_Options_ISD)         :: opt_sdc_isd

  namelist/time_integration_prm/ opt_pre_euler, opt_pre_isd &
                               , opt_sdc_euler, opt_sdc_euler_bepar &
                               , opt_sdc_isd

  real(RNP) :: c_min =   0  ! min CFL number
  real(RNP) :: c_max =  10  ! max CFL number
  real(RNP) :: d_min =  -5  ! min diffusion number
  real(RNP) :: d_max =  10  ! max diffusion number
  integer   :: nc    =  80  ! num intervals along imaginary axis
  integer   :: nd    =  80  ! num intervals along real axis

  namelist /input/ c_min, c_max, d_min, d_max, nc, nd

  real(RNP)   , allocatable :: c(:) ! convective CFL numbers (0:nc)
  real(RNP)   , allocatable :: d(:) ! diffusive  CFL numbers (0:nd)

  complex(RNP), allocatable :: lambda(:,:)    ! λ = -d + ic = λᵣ + iλᵢ
  real(RNP)   , allocatable :: a(:,:,:)       ! amplification: a = |u|
  real(RNP)   , allocatable :: e(:,:,:)       ! error        : e = |u -u_ex|
  real(RNP)                 :: c0, delta_c, delta_d, tol

  ! MLSDC
  type(DQ_MLSDC_1D)          :: mlsdc
  type(DQ_MLSDC_Options_1D)  :: mlsdc_opt

  type(DQ_MLSDC_Variable_1D) :: u_h, u_x
  real(RNP)                  :: t_0, t_1, dt
  real(RNP)                  :: err_max, r_max
  logical                    :: exists
  integer                    :: io, stat
  integer                    :: l, i, j, k
  character(len=80)          :: filename
  real(RNP)                  :: dt_f, int_f, int_c

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'MLSDC Component Testprogram Dahlquist'

  ! identify case
  call get_command_argument(1, case_name, status=stat)
  if (stat /= 0 .or. len_trim(case_name) == 0) then
    call get_command_argument(0, case_name, status=stat)
  end if

  case_file = trim(case_name) // '.prm'
  inquire(file=case_file, exist=exists)
  if (exists) then
    open(newunit=io, file=case_file)
    read(io, nml = discretization_prm)
    read(io, nml = time_integration_prm)
    read(io, nml = input)
    close(io)
  end if

  ! initialize parameters ......................................................

  allocate( c(0:nc), d(0:nd), lambda(0:nd,0:nc)         &
          , a(0:nd,0:nc,n_level), e(0:nd,0:nc,n_level))

  ! cfl numbers
  delta_c = (c_max - c_min) / max(nc,1)
  do j = 0, nc
    c(j) = c_min + j * delta_c
  end do

  ! diffusion numbers
  delta_d = (d_max - d_min) / max(nd,1)
  do i = 0, nd
    d(i) = d_min + i * delta_d
  end do

  ! λ = -d + ic
  do j = 0, nc
  do i = 0, nd
    lambda(i,j) =  -d(nd - i) + (ZERO,ONE) * c(j)
  end do
  end do

  ! MLSDC options
  mlsdc_opt = DQ_MLSDC_Options_1D(n_level)
  do l = 1, n_level
    mlsdc_opt % p_time  (l)   = p_time(l)
    mlsdc_opt % n_time  (l)   = n_time(l)
    mlsdc_opt % n_stage_p (l) = n_stage_p(l)
    mlsdc_opt % n_stage_c (l) = n_stage_c(l)
  end do
  write(*,*)
  write(*,'(A,99I5)') 'n_cycle  = ', n_cycle
  write(*,'(A,99I5)') 'n_coarse = ', n_coarse
  write(*,'(A,99I5)') 'n_level  = ', mlsdc_opt % n_level
  write(*,'(A,99I5)') 'p_time   = ', mlsdc_opt % p_time
  write(*,'(A,99I5)') 'n_time   = ', mlsdc_opt % n_time
  write(*,*)

  ! predictor options (ISD1, so far)
  select case(time_method)
  case(1)
    allocate(DQ_TimeIntegrator_Options_Euler :: opt_pre)
    opt_pre  =  opt_pre_euler
  case(3)
    allocate(DQ_TimeIntegrator_Options_ISD :: opt_pre)
    opt_pre  =  opt_pre_isd
  end select

  ! SDC options (ISD1 and Euler so far)
  select case(sdc_method)
  case(1)
    allocate(DQ_MLSDC_Corrector_Options_Euler :: opt_sdc)
    opt_sdc = opt_sdc_euler
  case(2)
    allocate(DQ_MLSDC_Corrector_Options_ISD :: opt_sdc)
    opt_sdc = opt_sdc_isd
  case(3)
    allocate(DQ_MLSDC_Corrector_Options_Euler_BEPar :: opt_sdc)
    opt_sdc = opt_sdc_euler_bepar
  end select

  ! MLSDC data structure
  mlsdc = DQ_MLSDC_1D(mlsdc_opt, opt_pre, opt_sdc)

  ! MLSDC variables
  u_h = DQ_MLSDC_Variable_1D(mlsdc)
  u_x = DQ_MLSDC_Variable_1D(mlsdc)

  ! testing of the components ..................................................

  t_0 = 0.0
  t_1 = t_0 + dt_slab

  ! test of the v-cycle ........................................................

  do j = 0, nc
  do i = 0, nd

  ! set initial condition
  do l = 1, n_level
    u_h%level(l)%val(0,1) = (1.0,0.0)
  end do

  select case(mg_start)
  case(0)

    ! Constant start
    do l = n_level, 1, -1
      associate( u_h => u_h % level(l) % val )

        ! set constant
        u_h(:,:) = (1.0,0.0)

      end associate
    end do

  case(1)

    ! Predictor start
    do l = n_level, 1, -1
      associate( u_h => u_h % level(l) % val )

        ! predictor on all levels
        call mlsdc % level(l) % ApplyPredictor(lambda(i,j), dt_slab, u_h)

      end associate
    end do

  case(2)

    ! Cascade start
    call DQ_MLSDC_Cascade_1D(mlsdc, incremental, lambda(i,j), dt_slab, 1, u_h)

  case(3)

    ! FMG start
    call DQ_MLSDC_FMG_1D( &
        mlsdc, incremental, lambda(i,j), dt_slab, 1, 1, n_coarse, u_h)

  end select

    ! enter v cycle
    call DQ_MLSDC_V_Cycle_1D( &
        mlsdc, incremental, lambda(i,j), dt_slab, 1, 1, n_coarse, n_cycle, u_h)

    do l = n_level, 1, -1
      associate( u_x    => u_x   % level(l) % val    &
               , u_h    => u_h   % level(l) % val    &
               , m_time => mlsdc % level(l) % m_time )

        a(i,j,l) = abs(u_h(m_time,1))
        u_x(0,1) = (1.0, 0.0)
        call GetExactSolution(mlsdc%level(l), t_0, t_1, lambda(i,j), u_x)
        e(i,j,l) = abs(u_h(m_time,1) - u_x(m_time,1))
      end associate
    end do

  end do
  end do

  do l = n_level, 1, -1
    err_max = maxval(e(:,:,l))
    write(*,'(2X,A,I3,A,ES10.3)') 'sdc error on level ',l,': err_max =',err_max

  ! save results ...............................................................

    ! x = Re(λ) = -d
    open(newunit=io, file='lambda_re.dat')
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
    write(filename, '(A,I0,A)') 'amplification_level_', l, '.dat'
    open(newunit=io, file = trim(filename))
    do j = 0, nc
      do i = 0, nd
        write(io, '(ES16.8E3,1X)', advance='NO') a(i,j,l)
      end do
      write(io,*)
    end do
    close(io)

    ! error
    write(filename, '(A,I0,A)') 'error_level_', l, '.dat'
    open(newunit=io, file = trim(filename))
    do j = 0, nc
      do i = 0, nd
        write(io, '(ES16.8E3,1X)', advance='NO') e(i,j,l)
      end do
      write(io,*)
    end do
    close(io)
  end do

contains

  !---------------------------------------------------------------------------
  !> Gets the exact solution of for a given mesh level

  subroutine GetExactSolution(level, t_0, t_1, lambda, u)
    class(DQ_MLSDC_Level_1D), intent(in)  :: level           !< space-time level
    real(RNP),                intent(in)  :: t_0             !< initial time
    real(RNP),                intent(in)  :: t_1             !< final time
    complex(RNP),             intent(in)  :: lambda
    complex(RNP),             intent(out) :: u(0:,:)

    real(RNP), allocatable :: t(:)
    real(RNP) :: dt
    integer   :: i, j, k, l, nd, nc, mt, nt

    mt = ubound(u,1)
    nt = ubound(u,2)
    dt = (t_1 - t_0) / nt

    allocate(t(0:mt))

    associate(dq_sdc => level % sdc)

      do l = 1, nt
        t(0:) = dq_sdc % SubintervalPoints(t_0 + (l-1)*dt, dt)
        do k = 0, mt
          u(k,l) = exp(lambda * (t(k)-t_0)) ! * u_0 = (1,0)
        end do
      end do

    end associate

  end subroutine GetExactSolution

  !=============================================================================

end program Dahlquist_MLSDC
