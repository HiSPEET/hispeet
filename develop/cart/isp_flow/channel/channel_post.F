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

!> summary:  Postprocessing for Channel flow case
!>
!> Use the first argument of the invoking command to specify the base name of
!> the control file. If omitted, the program looks for `channel_post.prm`.
!===============================================================================

!!! suggestion for parallelization
!!! - move data of xy planes to different MPI processes
!!! – perform planform averaging and FFTs in parallel, i.e. individually
!!!   for each MPI process holding one or more xy planes
!!! - collect reduced data on MPI master
!!! – perform final computations and output with MPI master

program Channel_Post

  use, intrinsic :: IEEE_Arithmetic, only: ieee_is_nan

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, TWO, PI, HALF
  use OpenMP_Binding
  use FFTW_Binding
  use Execution_Control
  use Array_Assignments
  use Array_Reductions
  use Gauss_Jacobi
  use Standard_Element_Operators__1D
  use XMPI

  use ISP_Flow_Problem__Test_Suite

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP .............................................................

  type(MPI_Comm) :: comm         ! MPI communicator
  integer        :: rank         ! local MPI rank
  integer        :: n_proc       ! number of MPI processes
  integer        :: n_thread     ! number of OpenMP threads

  type(MPI_Status) :: status

  ! control ....................................................................

  character(len=80) :: control_file

  character(len=80) :: flow_case    = 'channel'
  character(len=80) :: solver_file  = 'isp_flow__sdc_test'

  logical :: use_symm = .false. ! calculate the average useing the symmetry of problem

  ! parameters read from control file
  namelist /control/ flow_case, solver_file,    &
                     use_symm

  ! discretization .............................................................

  ! parameters read from solver file
  integer   :: np(3)        ! number of partitions in directions 1:3
  integer   :: ep(3)        ! elements per partition and direction
  integer   :: po_u         ! polynomial order of u
  integer   :: po_p         ! -- ignored --
  integer   :: po_q         ! -- ignored --
  real(RNP) :: penalty      ! -- ignored --
  logical   :: adjust_dx    ! -- ignored --

  namelist /discretization/ np, ep, po_u, po_p, po_q, penalty, adjust_dx

  ! derived parameters
  real(RNP) :: dx(3)        ! spacing in directions 1:3
  integer   :: ip(3)        ! 3D partition ID, 0 ≤ ip(k) < np(k), k = 1:3
  integer   :: ne           ! number of elements in the local partition
  integer   :: ni           ! interpolated points per element and direction
  integer   :: no           ! originally points per element and direction
  integer   :: nl(3)        ! local points per direction
  integer   :: ng(3)        ! global points per direction

  real(RNP) :: l, w, h      ! length, width and hight of channel

  ! flow problem ...............................................................

  class(FlowProblem), allocatable :: problem
  real(RNP) :: t

  ! operators ..................................................................

  type(StandardElementOperators_1D), save :: standard_op
  type(C_PTR) :: plan

  ! variables ..................................................................

  ! computed solution from input flow data
  real(RNP), allocatable, save :: uc(:,:,:,:,:)       ! computed solution
  real(RNP), allocatable, save :: q_avg(:,:,:,:,:)    ! time-averaged quantities
  real(RNP), allocatable, save :: x3(:,:)             ! z-coordinates [0,h]

  ! interpolated solution
  real(RNP), allocatable, save :: ui(:,:,:,:,:)       ! interpolated quantities
  real(RNP), allocatable, save :: q_avg_i(:,:,:,:,:)  ! interpolated time-averaged quantities
  real(RNP), allocatable, save :: xi(:)               ! interpolation points in [-1,1]

  ! global/transformed quantities on regular Cartesian mesh (rank 0 only)
  real    (C_DOUBLE)        , allocatable, save :: ug(:,:,:,:,:)
  real    (C_DOUBLE)        , allocatable, save :: q_avg_g(:,:,:,:,:)
  complex (C_DOUBLE_COMPLEX), allocatable, save :: vt_1d_1(:,:,:,:,:)
  complex (C_DOUBLE_COMPLEX), allocatable, save :: vt_1d_2(:,:,:,:,:)
  complex (C_DOUBLE_COMPLEX), allocatable, save :: vt_2d(:,:,:,:,:)

  ! 1D and 2D premultiplied energy spectra
  real    (C_DOUBLE)        , allocatable, save :: e_1d_1(:,:,:,:)
  real    (C_DOUBLE)        , allocatable, save :: e_1d_2(:,:,:,:)
  real    (C_DOUBLE)        , allocatable, save :: e_2d(:,:,:,:,:)

  ! plane-averaged velocity and derivatives
  real(RNP), allocatable, save :: v_mean(:,:,:)       ! time- and plane-averaged velocity
  real(RNP), allocatable, save :: dv_mean_dz(:,:,:)   ! time- and plane-averaged velocity derivatives

  ! Computed friction Reynolds number
  real(RNP) :: Re_tau_comp

  ! Reynolds stress tensor (rank 0 only)
  real(RNP), allocatable, save :: r(:,:,:,:,:)        ! Reynolds stress tensor on global field
  real(RNP), allocatable, save :: r_planeavg(:,:,:)   ! plane-averaged Reynolds stress tensor

  logical   :: exists
  integer   :: prm, stat
  integer   :: eg3                  ! number of elements in z-direction
  integer   :: e, e3
  integer   :: c
  integer   :: i, j, k
  integer   :: k1, k2

  !-----------------------------------------------------------------------------
  ! initialization

  ! MPI ........................................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  ! control and discretization parameters ......................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Postprocessing for Channel problem'
    write(*,*)

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'channel_post'
    end if

    ! read control file
    inquire(file=trim(control_file)//'.prm', exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(control_file)//'.prm'
      open(newunit=prm, file=trim(control_file)//'.prm')
      read(prm, nml=control)
      close(prm)
    end if

    ! read solver file
    inquire(file=trim(solver_file)//'.prm', exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(solver_file)//'.prm'
      open(newunit=prm, file=trim(solver_file)//'.prm')
      read(prm, nml=discretization)
      close(prm)
    end if

    write(*,'(2X,A,1X,A)')  'flow case   = ', trim(flow_case)
    write(*,'(2X,A,1X,A)')  'solver file = ', trim(solver_file)
    write(*,'(2X,A,1X,I0)') 'n_proc      = ', n_proc
    write(*,'(2X,A,1X,I0)') 'n_thread    = ', n_thread
    write(*,*)

  end if

  ! globalize parameters
  call XMPI_Bcast(np  , 0, comm)
  call XMPI_Bcast(ep  , 0, comm)
  call XMPI_Bcast(po_u, 0, comm)

  ! flow problem ...............................................................

  allocate(FlowProblem_Channel :: problem)

  call problem % SetProblem(flow_case, comm)

  ! derived parameters .........................................................

  l     = problem % x1(1) - problem % x0(1)
  w     = problem % x1(2) - problem % x0(2)
  h     = problem % x1(3) - problem % x0(3)
  dx    = (problem % x1 - problem % x0) / (np * ep)
  ne    = product(ep)
  no    = po_u + 1
  ni    = po_u + 1
  nl(1) = ni        * ep(1)
  nl(2) = ni        * ep(2)
  nl(3) = no * ep(3)
  ng    = nl * np
  eg3   = ep(3)*np(3)

  ! 3D partition ID
  ip = CartPID(rank, np)

  if (rank == 0) then
    write(*,*)
    write(*,'(A,3I5)') 'np = ', np
    write(*,'(A,3I5)') 'ip = ', ip
    write(*,'(A,3I5)') 'ep = ', ep
    write(*,'(A,3I5)') 'ni = ', ni
    write(*,'(A,3I5)') 'nl = ', nl
    write(*,'(A,3I5)') 'ng = ', ng
  end if

  ! variables and operators ....................................................

  standard_op = StandardElementOperators_1D(po_u)

  ! interpolation points of homogeneous directions
  allocate(xi(ni))
  do i = 1, ni
    xi(i) = (TWO * i - ONE)/ni - ONE
  end do

  ! global coordinates in heterogeneous direction
  if (rank == 0) then
    allocate( x3(0:po_u, eg3)  )
    do e3 = 1, eg3
    do k  = 0, po_u
      x3(k,e3) =  dx(3) * (e3 - HALF + HALF*standard_op%x(k))
    end do
    end do
  end if

  !-----------------------------------------------------------------------------
  ! Read and interpolate flow data

  ! load
  allocate(    uc(0:po_u, 0:po_u, 0:po_u, ne, problem%nc),  &
            q_avg(0:po_u, 0:po_u, 0:po_u, ne, 9)            )

  call LoadFlowVariables(flow_case, t, uc, q_avg)

  ! interpolate in x1 and x2 to equidistant nodes
  allocate(     ui(1:ni, 1:ni, 0:po_u, ne, problem%nc), &
           q_avg_i(1:ni, 1:ni, 0:po_u, ne, 9)  )

  call Interpolation2D(no, ni, ne, problem%nc, xi, uc, ui)
  call Interpolation2D(no, ni, ne, 9, xi, q_avg, q_avg_i)

  deallocate(uc, q_avg)

  !-----------------------------------------------------------------------------
  ! Transfer to data to rank 0

  ! Transfer velocity
  if (rank == 0) then
    allocate(ug(1:ng(1), 1:ng(2), 0:po_u, 1:eg3, problem%nc))
  end if

  if (rank > 0) then
    call TransferDataToRank0(status, np, nl, ep, ip, ui)
  else
    call TransferDataToRank0(status, np, nl, ep, ip, ui, ug)
  end if

  deallocate(ui)

  ! Transfer time-averaged quantities
  if (rank == 0) then
    allocate(q_avg_g(1:ng(1), 1:ng(2), 0:po_u,1:eg3 ,9))
  end if

  if (rank > 0) then
    call TransferDataToRank0(status, np, nl, ep, ip, q_avg_i)
  else
    call TransferDataToRank0(status, np, nl, ep, ip, q_avg_i, q_avg_g)
  end if

  deallocate(q_avg_i)

  !-----------------------------------------------------------------------------
  ! Fourier transformation of velocity

  if (rank == 0) then

    ! 1D .......................................................................

    !!! add description what happens here !!!

    ! streamwise
    allocate(vt_1d_1(0:ng(1)/2, 1:ng(2), 0:po_u, 1:eg3 , 3))
    plan = fftw_plan_dft_r2c_1d( ng(1), ug(:,1,0,1,1), vt_1d_1(:,1,0,1,1) &
                               , FFTW_ESTIMATE )
    do c  = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
    do j  = 1, ng(2)
        call fftw_execute_dft_r2c(plan, ug(:,j,k,e3,c), vt_1d_1(:,j,k,e3,c))
    end do
    end do
    end do
    end do

    call fftw_destroy_plan(plan)

    vt_1d_1 = vt_1d_1 / ng(1)

    ! spanwise
    allocate(vt_1d_2(1:ng(1), 0:ng(2)/2, 0:po_u, 1:eg3 , 3))
    plan = fftw_plan_dft_r2c_1d( ng(2), ug(1,:,0,1,1), vt_1d_2(1,:,0,1,1) &
                               , FFTW_ESTIMATE )

    do c  = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
    do i  = 1, ng(1)
        call fftw_execute_dft_r2c(plan, ug(i,:,k,e3,c), vt_1d_2(i,:,k,e3,c))
    end do
    end do
    end do
    end do

    call fftw_destroy_plan(plan)

    vt_1d_2 = vt_1d_2 / ng(2)

    ! 2D FFT ...................................................................

    !!! add description what happens here !!!

    allocate(vt_2d(0:ng(1)/2, 0:ng(2)-1, 0:po_u, 1:eg3 , 3))

    !!! would expect rather something like fftw_plan_dft_r2c_2d !!!

    plan = fftw_plan_dft_1d( ng(2)              &
                           , vt_1d_1(1,:,0,1,1) &
                           , vt_2d(1,:,0,1,1)   &
                           , -1                 &
                           , FFTW_ESTIMATE      )

    do c  = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
    do i  = 0, ng(1)/2
      call fftw_execute_dft(plan, vt_1d_1(i,:,k,e3,c), vt_2d(i,:,k,e3,c))
    end do
    end do
    end do
    end do

    call fftw_destroy_plan(plan)

    vt_2d = vt_2d / product(ng(1:2))

  end if

  !-----------------------------------------------------------------------------
  ! Mean velocities

  if (rank == 0) then

    allocate(v_mean(0:po_u,1:eg3,2))
    call PlaneAverage(ng(1), ng(2), no, eg3, 2, q_avg_g(:,:,:,:,1:2), v_mean)

  end if

  !-----------------------------------------------------------------------------
  ! derivatives of mean velocities

  if (rank == 0) then

    allocate(dv_mean_dz(0:po_u,1:eg3,2))
    call GetDerivative(no, eg3, standard_op%D, dx(3), v_mean(:,:,1), dv_mean_dz(:,:,1))
    call GetDerivative(no, eg3, standard_op%D, dx(3), v_mean(:,:,2), dv_mean_dz(:,:,2))

  end if

  !-----------------------------------------------------------------------------
  ! Computed friction Reynolds number

  if (rank == 0) then

    Re_tau_comp = (dv_mean_dz(0,1,1) - dv_mean_dz(no-1, eg3,1))/2

  end if

  !-----------------------------------------------------------------------------
  ! Reynolds stress

  if (rank == 0) then

    allocate( r(1:ng(1), 1:ng(2), 0:po_u, 1:eg3, 6) )

    do e3 = 1, eg3
      do k  = 0, po_u
      do j  = 1, ng(2)
      do i  = 1, ng(1)

        ! <u'u'>
        r(i,j,k,e3,1) = q_avg_g(i,j,k,e3,4) &
                      - q_avg_g(i,j,k,e3,1) * q_avg_g(i,j,k,e3,1)
        ! <u'v'>
        r(i,j,k,e3,2) = q_avg_g(i,j,k,e3,5) &
                      - q_avg_g(i,j,k,e3,2) * q_avg_g(i,j,k,e3,1)
        ! <u'w'>
        r(i,j,k,e3,3) = q_avg_g(i,j,k,e3,6) &
                      - q_avg_g(i,j,k,e3,3) * q_avg_g(i,j,k,e3,1)
        ! <v'v'>
        r(i,j,k,e3,4) = q_avg_g(i,j,k,e3,7) &
                      - q_avg_g(i,j,k,e3,2) * q_avg_g(i,j,k,e3,2)
        ! <v'w'>
        r(i,j,k,e3,5) = q_avg_g(i,j,k,e3,8) &
                      - q_avg_g(i,j,k,e3,3) * q_avg_g(i,j,k,e3,2)
        ! <w'w'>
        r(i,j,k,e3,6) = q_avg_g(i,j,k,e3,9) &
                      - q_avg_g(i,j,k,e3,3) * q_avg_g(i,j,k,e3,3)

      end do
      end do
      end do
    end do

    ! plane-average (arithmetic mean)
    allocate( r_planeavg(0:po_u, 1:eg3, 6) )

    call PlaneAverage(ng(1), ng(2), no, eg3, 6, r, r_planeavg)

  end if

  !-----------------------------------------------------------------------------
  ! energy spectra 1D and 2D

  if (rank == 0) then

  ! 1D Energy spectra ..........................................................

  ! streamwise
  allocate(e_1d_1(0:ubound(vt_1d_1,1), 0:po_u, 1:eg3, 1:3))

  e_1d_1 = 0

  do c = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
    do j  = 1, ng(2)
      do k1 = 0, ubound(vt_1d_1,1)
          e_1d_1(k1,k,e3,c) = e_1d_1(k1,k,e3,c) &
                                + HALF * abs(vt_1d_1(k1,j,k,e3,c))**2
      end do
    end do
    end do
    end do
  end do

  e_1d_1 = e_1d_1 / ng(2)

  ! spanwise
  allocate(e_1d_2(0:ubound(vt_1d_2,2), 0:po_u, 1:eg3, 1:3))

  e_1d_2 = 0

  do c = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
      do k2 = 0, ubound(vt_1d_2,2)
        do i  = 1, ng(1)
          e_1d_2(k2,k,e3,c) = e_1d_2(k2,k,e3,c) &
                                + HALF * abs(vt_1d_2(i,k2,k,e3,c))**2
        end do
      end do
    end do
    end do
  end do

  e_1d_2 = e_1d_2 / ng(1)

  ! 2D Energy spectra ..........................................................

  allocate(e_2d(0:ubound(vt_2d,1), 0:ubound(vt_2d,2), 0:po_u, 1:eg3, 1:3))

  do c = 1, 3
    do e3 = 1, eg3
    do k  = 0, po_u
      do k2 = 0, ubound(vt_2d,2)
      do k1 = 0, ubound(vt_2d,1)
        e_2d(k1,k2,k,e3,c) = HALF * abs(vt_2d(k1,k2,k,e3,c))**2
      end do
      end do
    end do
    end do
  end do

  end if

  !-----------------------------------------------------------------------------
  ! Apply symmetry with respect to midplane

  if (use_symm) then
    if (rank == 0) then

    call SymmetryAverage(no, eg3,  1, v_mean(:,:,1))
    call SymmetryAverage(no, eg3, -1, v_mean(:,:,2))
    call SymmetryAverage(no, eg3, -1, dv_mean_dz(:,:,1))
    call SymmetryAverage(no, eg3,  1, dv_mean_dz(:,:,2))
    call SymmetryAverage(no, eg3,  1, r_planeavg(:,:,1))
    call SymmetryAverage(no, eg3, -1, r_planeavg(:,:,2))
    call SymmetryAverage(no, eg3, -1, r_planeavg(:,:,3))
    call SymmetryAverage(no, eg3,  1, r_planeavg(:,:,4))
    call SymmetryAverage(no, eg3, -1, r_planeavg(:,:,5))
    call SymmetryAverage(no, eg3,  1, r_planeavg(:,:,6))
    do k1 = 0, ubound(vt_1d_1,1)
      call SymmetryAverage(no, eg3,  1, e_1d_1(k1,:,:,1))
      call SymmetryAverage(no, eg3,  1, e_1d_1(k1,:,:,2))
      call SymmetryAverage(no, eg3,  1, e_1d_1(k1,:,:,3))
    end do
    do k2 = 0, ubound(vt_1d_2,2)
      call SymmetryAverage(no, eg3,  1, e_1d_1(k2,:,:,1))
      call SymmetryAverage(no, eg3,  1, e_1d_1(k2,:,:,2))
      call SymmetryAverage(no, eg3,  1, e_1d_1(k2,:,:,3))
    end do
    end if
  end if

  !-----------------------------------------------------------------------------
  ! Evaluation

  if (rank == 0) then

    ! Computed friction Reynolds number
    write(*,'(/,A,F7.2)') 'Re_tau = ', Re_tau_comp


    ! Mean velocities their derivatives
    write(*,'(/,A,/)') 'Means'
    write(*,'(A, 1X)',advance='NO') '#'
    write(*,'(A, 3X)',advance='NO') '    z     '
    write(*,'(A, 3X)',advance='NO') '    z+    '
    write(*,'(A, 3X)',advance='NO') '   Umean  '
    write(*,'(A, 3X)',advance='NO') '   Vmean  '
    write(*,'(A, 3X)',advance='NO') ' dUmean/dz'
    write(*,'(A, 3X)',advance='NO') ' dVmean/dz'
    write(*,*)
    do e3 = 1, eg3
      do k = 0, po_u
        write(*,'(ES12.5,1X)',advance='NO') x3(k,e3)
        write(*,'(ES12.5,1X)',advance='NO') x3(k,e3)*Re_tau_comp
        write(*,'(2(ES12.5,1X))',advance='NO') v_mean(k,e3,:)
        write(*,'(2(ES12.5,1X))',advance='NO') dv_mean_dz(k,e3,:)
        write(*,*)
      end do
    end do

    ! Reynolds Stress
    write(*,'(/,A,/)') 'Reynolds Stress'
    write(*,'(A, 1X)',advance='NO') '#'
    write(*,'(A, 3X)',advance='NO') '    z     '
    write(*,'(A, 3X)',advance='NO') '    z+    '
    write(*,'(A, 3X)',advance='NO') '   R_uu   '
    write(*,'(A, 3X)',advance='NO') '   R_uv   '
    write(*,'(A, 3X)',advance='NO') '   R_uw   '
    write(*,'(A, 3X)',advance='NO') '   R_vv   '
    write(*,'(A, 3X)',advance='NO') '   R_vw   '
    write(*,'(A, 3X)',advance='NO') '   R_ww   '
    write(*,*)
    do e3 = 1, eg3
      do k  = 0, po_u
        write(*,'(ES12.5,1X)',advance='NO') x3(k,e3)
        write(*,'(ES12.5,1X)',advance='NO') x3(k,e3)*Re_tau_comp
        write(*,'(6(ES12.5,1X))',advance='NO') r_planeavg(k,e3,:)
        write(*,*)
      end do
    end do

    ! Energy spectra
    write(*,'(/,A,/)') 'Energy spectra in x-direction'
    do e3 = 1, eg3
    do k  = 0, po_u
      write(*,'(/,A,ES12.5,/)') 'z+ = ', x3(k,e3)*Re_tau_comp
      write(*,'(A, 1X)',advance='NO') '#'
      write(*,'(A, 3X)',advance='NO') '  k  '
      write(*,'(A, 3X)',advance='NO') '   E_xx   '
      write(*,'(A, 3X)',advance='NO') '   E_xy   '
      write(*,'(A, 3X)',advance='NO') '   E_xz   '
      write(*,*)
      do k1 = 0, ubound(vt_1d_1,1)
          if (sum(e_1d_1(k1,k,e3,:)) < epsilon(1.0)) cycle
          write(*,'(F6.1,1X)',advance='NO') k1*2*PI/l
          write(*,'(3(ES12.5,1X))',advance='NO') e_1d_1(k1,k,e3,:)
          write(*,*)
      end do
    end do
    end do

    write(*,'(/,A,/)') 'Energy spectra in y-direction'
    do e3 = 1, eg3
    do k  = 0, po_u
      write(*,'(/,A,ES12.5,/)') 'z+ = ', x3(k,e3)*Re_tau_comp
      write(*,'(A, 1X)',advance='NO') '#'
      write(*,'(A, 3X)',advance='NO') '  k  '
      write(*,'(A, 3X)',advance='NO') '   E_yx   '
      write(*,'(A, 3X)',advance='NO') '   E_yy   '
      write(*,'(A, 3X)',advance='NO') '   E_yz   '
      write(*,*)
      do k2 = 0, ubound(vt_1d_2,2)
          if (sum(e_1d_2(k2,k,e3,:)) < epsilon(1.0)) cycle
          write(*,'(F6.1,1X)',advance='NO') k2*2*PI/w
          write(*,'(3(ES12.5,1X))',advance='NO') e_1d_2(k2,k,e3,:)
          write(*,*)
      end do
    end do
    end do

  end if

  !---------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!-----------------------------------------------------------------------------
!> Map given MPI rank to 3D Cartesian partition ID

function CartPID(rk, np) result(ip)
  integer, intent(in)  :: rk     !< MPI rank
  integer, intent(in)  :: np(3)  !< number of partitions in directions 1:3
  integer              :: ip(3)  !< triple partition ID, 0 ≤ ip < np

  ip(3) = rk / (np(1) * np(2))
  ip(2) = (rk - np(1) * np(2) * ip(3)) / np(1)
  ip(1) = rk - np(1) * (ip(2) + np(2) * ip(3))

end function CartPID


!-----------------------------------------------------------------------------
!> Transfer local data from all ranks to global data field on rank 0

subroutine TransferDataToRank0(status, np, nl, ep, ip, ul, ug)
  type(MPI_Status),    intent(in)    :: status
  integer,             intent(in)    :: np(3)           !< number of partitions per direction
  integer,             intent(in)    :: nl(3)           !< local points per direction
  integer,             intent(in)    :: ep(3)           !< elements per partition and direction
  integer,             intent(in)    :: ip(3)           !< 3D partition ID, 0 ≤ ip(k) < np(k), k = 1:3
  real(RNP),           intent(in)    :: ul(:,:,:,:,:)   !< local data field
  real(RNP), optional, intent(inout) :: ug(:,:,:,:,:)   !< global data field

  integer :: ip_i(3)       ! 3D partition ID of ith partition
  integer :: i

  if (rank > 0) then

    call MPI_Send(ul, size(ul), MPI_REAL_RNP, 0, 0, MPI_COMM_WORLD)

  else

    ! own contribution
    call CopyLocalToGlobal(nl, ep, ip, ul, ug)

    do i = 1, product(np)-1

      ! receive data from partition i
      call MPI_Recv(ul, size(ul), MPI_REAL_RNP, i, 0, MPI_COMM_WORLD, status)

      ! copy into global array
      ip_i = CartPID(i, np)
      call CopyLocalToGlobal(nl, ep, ip_i, ul, ug)

    end do

  end if

end subroutine TransferDataToRank0

!-----------------------------------------------------------------------------
!> Copy local data to global cartesian array

subroutine CopyLocalToGlobal(nl, ep, ip, ul, ug)
  integer,   intent(in)    :: nl(3)           !< local points per direction
  integer,   intent(in)    :: ep(3)           !< elements per partition and direction
  integer,   intent(in)    :: ip(3)           !< 3D partition ID, 0 ≤ ip(k) < np(k), k = 1:3
  real(RNP), intent(in)    :: ul(:,:,:,:,:)   !< local data field
  real(RNP), intent(inout) :: ug(:,:,:,:,:)   !< global data field

  integer :: ni
  integer :: nc
  integer :: e1, e2, e3, eg3
  integer :: i0, il, ig
  integer :: j0, jl, jg
  integer :: kl
  integer :: c

  ni = size(ul,1)
  nc = size(ul,5)

  ! global offsets
  i0 = nl(1) * ip(1)
  j0 = nl(2) * ip(2)
  do c = 1, nc
    do e3 = 1, ep(3)
    do e2 = 1, ep(2)
    do e1 = 1, ep(1)

      ! element ID
      e  = e1 + ep(1) * (e2 - 1 + ep(2) * (e3 - 1))

      ! global element ID in z-direction
      eg3 = e3 + (ip(3)*ep(3))

      do kl = 1, po_u + 1
      do jl = 1, ni
      do il = 1, ni

        ! global point indices
        ig = i0 +  ni      * (e1 - 1) + il
        jg = j0 +  ni      * (e2 - 1) + jl

        ug(ig,jg,kl,eg3,c) = ul(il,jl,kl,e,c)

      end do
      end do
      end do

    end do
    end do
    end do
  end do

end subroutine CopyLocalToGlobal

!-----------------------------------------------------------------------------
!> Loads the flow variables from file

subroutine LoadFlowVariables(flow_case, t, u, q_avg)
  character(len=*),    intent(in)  :: flow_case        !< flow case
  real(RNP),           intent(out) :: t                !< time
  real(RNP),           intent(out) :: u(:,:,:,:,:)     !< flow variables
  real(RNP), optional, intent(out) :: q_avg(:,:,:,:,:) !< time-averaged quantities

  integer ::  n_avg                !< current number of time-averages

  integer :: unit, ios
  character(len=100) :: file

  !$omp single

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'

  open(newunit=unit, file=file, status='OLD', form='UNFORMATTED', iostat=ios)
  if (ios == 0) then
    read(unit) t
    read(unit) u
    read(unit) n_avg
    if (present(q_avg)) then
      read(unit) q_avg
    else
      call SetArray(q_avg, ZERO, multi = .true.)
    end if
    close(unit)
  else
    call Error('LoadFlowVariables', 'found no matching input file')
  end if

  !$omp end single

end subroutine LoadFlowVariables

!-----------------------------------------------------------------------------
!> Interpolation of 2D mesh variables

subroutine Interpolation2D(no, ni, ne, nc, xi, uo, ui)
  integer,   intent(in)    :: no                             !< number of original mesh points per element
  integer,   intent(in)    :: ni                             !< number of interpolated mesh points per element
  integer,   intent(in)    :: ne                             !< number of elements per patition
  integer,   intent(in)    :: nc                             !< number of components (flow variables)
  real(RNP), intent(in)    :: xi(0:ni-1)                     !< interpolation points in [-1,1]
  real(RNP), intent(in)    :: uo(0:no-1,0:no-1,0:no-1,ne,nc) !< original mesh variable
  real(RNP), intent(inout) :: ui(0:ni-1,0:ni-1,0:no-1,ne,nc) !< interpolated mesh variable

  !---------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: At(:,:) !< transposed 1D interpolation operator
  real(RNP), allocatable :: z2(:,:) !< z2 = AxI uo
  real(RNP) :: tmp
  integer :: i, j, k, p, po, pi

  associate(xo => standard_op%x)

    po  = ubound(xo,1)
    pi  = ubound(xi,1)

    allocate( At(0:po, 0:pi), &
              z2(0:po, 0:pi)  )

    do j = 0, po
    do i = 0, pi
      At(j,i) = LobattoPolynomial(j, xo, xi(i))
    end do
    end do

    do c = 1, nc
    do e = 1, ne
    do k = 0, no-1

      ! z2 = AxI uo ........................................................

      do j = 0, pi
        do i = 0, po
          tmp = 0
            do p = 0, po
              tmp = tmp + At(p,j) * uo(i,p,k,e,c)
            end do
          z2(i,j) = tmp
        end do
      end do

      ! ui = IxA z2 ........................................................

      do j = 0, pi
        do i = 0, pi
          tmp = 0
          do p = 0, po
            tmp = tmp + At(p,i) * z2(p,j)
          end do
          ui(i,j,k,e,c) = tmp
        end do
      end do

    end do
    end do
    end do

    deallocate(At, z2)

  end associate

end subroutine Interpolation2D

!-----------------------------------------------------------------------------
!> Perform a plane-average of equidistant nodes (arithmetic mean)

subroutine PlaneAverage(ng1, ng2, n3, eg3, nc, ug, u_planeavg)
  integer  , intent(in)    :: ng1                       !< number of global points in x1-direction
  integer  , intent(in)    :: ng2                       !< number of global points in x2-direction
  integer  , intent(in)    :: n3                        !< number of local points in x3-direction
  integer  , intent(in)    :: eg3                       !< number of elements in x3-direction
  integer  , intent(in)    :: nc                        !< number of components (flow variables)
  real(RNP), intent(in)    :: ug(ng1,ng2,0:n3-1,eg3,nc) !< original mesh variables
  real(RNP), intent(inout) :: u_planeavg(0:n3-1,eg3,nc) !< plane-averaged mesh variables

  integer :: i, j, k, e3, c

  if (rank == 0) then

    u_planeavg = 0

      do c  = 1, nc
      do e3 = 1, eg3
      do k  = 0, n3-1

        do j  = 1, ng2
        do i  = 1, ng1

          u_planeavg(k,e3,c) = u_planeavg(k,e3,c) + ug(i,j,k,e3,c)

        end do
        end do

      end do
      end do
      end do

      u_planeavg = u_planeavg / (ng1*ng2)

  end if

end subroutine PlaneAverage

!-----------------------------------------------------------------------------
!> Average of upper and lower channel half

subroutine SymmetryAverage(n3, eg3, s, u)
  integer  , intent(in)    :: n3                   !< number of local points in x3-direction
  integer  , intent(in)    :: eg3                  !< number of elements in x3-direction
  integer  , intent(in)    :: s                    !< sign (+1 for axial and -1 central symmetry)
  real(RNP), intent(inout) :: u(0:n3-1,eg3)        !< mesh variables

  real(RNP) :: u_sym(0:n3-1,eg3)
  integer :: k, e3

  if (rank == 0) then

    do e3 = 1, eg3
    do k  = 0, n3-1

      u_sym(k,e3) = (u(k,e3) + s * u(n3-1-k, eg3-e3+1)) / 2

    end do
    end do

    u = u_sym

  end if

end subroutine SymmetryAverage

!-----------------------------------------------------------------------------
!> Get derivative of plane-averaged variable in z direction

subroutine GetDerivative(n3, eg3, Ds, dx, u, du_dz)
  integer  , intent(in)    :: n3            !< number of local points in x3-direction
  integer  , intent(in)    :: eg3           !< number of elements in x3-direction
  real(RNP), intent(in)    :: Ds(n3,n3)     !< 1D standard diff matrix
  real(RNP), intent(in)    :: dx            !< element extension in z
  real(RNP), intent(in)    :: u(n3,eg3)     !< original mesh variables
  real(RNP), intent(inout) :: du_dz(n3,eg3) !< derivatived mesh variables

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP) :: A(n3,n3)
  real(RNP) :: g, tmp
  integer   :: p, k, e3

  !!! using A is more confusing tha helpful here, its only 1D ;) !!!

  A = transpose(Ds)

  g = 2 / dx

  do e3 = 1, eg3
    do k = 1, n3

      tmp = 0

      do p = 1, n3
        tmp = tmp + A(p,k) * u(p,e3)
      end do

      du_dz(k,e3) = g * tmp

    end do
  end do

end subroutine GetDerivative

!=============================================================================

end program Channel_Post
