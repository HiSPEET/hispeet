!> summary:  Postprocessing for Transition_TG flow case
!> author:   Joerg Stiller
!> date:     2020/10/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Use the first argument of the invoking command to specify the base name of
!> the control file. If omitted, the program looks for `transition_tg_post.prm`.
!===============================================================================

program Transition_TG_Post

  use, intrinsic :: IEEE_Arithmetic, only: ieee_is_nan

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, TWO, PI
  use OpenMP_Binding
  use FFTW_Binding
  use Execution_Control
  use Array_Assignments
  use Array_Reductions
  use Standard_Operators_1D
  use Embedded_Interpolation_3D
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

  ! control ....................................................................

  character(len=80) :: control_file
  character(len=80) :: flow_case    = 'transition_tg'
  character(len=80) :: solver_file  = 'isp_flow__sdc_test'

  integer :: kc(3) = -1   ! cosine wave number for testing (< 0 to skip)
  integer :: ks(3) = -1   ! sine   wave number for testing (< 0 to skip)

  ! parameters read from control file
  namelist /control/ flow_case, solver_file, kc, ks

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
  integer   :: nl(3)        ! local points per direction
  integer   :: ng(3)        ! global points per direction

  ! flow problem ...............................................................

  class(FlowProblem), allocatable :: problem
  real(RNP) :: t

  ! operators ..................................................................

  type(StandardOperators1D)    , save :: standard_op
  type(InterpolationOperator3D), save :: interpol_op

  type(C_PTR) :: plan

  ! variables ..................................................................

  real(RNP), allocatable, save :: uc(:,:,:,:,:) ! computed solution
  real(RNP), allocatable, save :: vi(:,:,:,:,:) ! interpolated velocity
  real(RNP), allocatable, save :: xi(:)         ! interpolation points in [-1,1]

  ! interpolated/transformed velocity on regular Cartesian mesh
  real    (C_DOUBLE)        , allocatable, save :: wi(:,:,:,:)
  complex (C_DOUBLE_COMPLEX), allocatable, save :: wt(:,:,:,:)

  logical   :: exists
  integer   :: prm, stat
  integer   :: e, e1, e2, e3
  integer   :: i, i0, i1, i2, ig
  integer   :: j, j0, j1, j2, jg
  integer   :: k, k0, k1, k2, kg

  complex(RNP) :: z

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
    write(*,'(A)') 'Postprocessing of TG transisition problem'
    write(*,*)

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'transition_tg_post'
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

  end if

  ! globalize parameters
  call XMPI_Bcast(np  , 0, comm)
  call XMPI_Bcast(ep  , 0, comm)
  call XMPI_Bcast(po_u, 0, comm)

  ! derived parameters
  dx = (problem % x1 - problem % x0) / (np * ep)
  ne = product(ep)
  ni = po_u + 1
  nl = ni * ep
  ng = nl * np

  ! 3D partition ID
  ip(3) = rank / (np(1) * np(2))
  ip(2) = (rank - np(1) * np(2) * ip(3)) / np(1)
  ip(1) = rank - np(1) * (ip(2) + np(2) * ip(3))

  if (rank == 0) then
    write(*,*)
    write(*,'(A,3I5)') 'np = ', np
    write(*,'(A,3I5)') 'ip = ', ip
    write(*,'(A,3I5)') 'ep = ', np
    write(*,'(A,3I5)') 'ni = ', ni
    write(*,'(A,3I5)') 'nl = ', nl
    write(*,'(A,3I5)') 'ng = ', ng
  end if

  ! flow problem ...............................................................

  allocate(FlowProblem_Transition_TG :: problem)
  call problem % SetProblem(flow_case, comm)

  ! variables and operators ....................................................

  allocate(uc(0:po_u, 0:po_u, 0:po_u, ne, problem%nc))
  call LoadFlowVariables(flow_case, t, uc)

  allocate(xi(ni))
  do i = 1, ni
    xi(i) = (TWO * i - ONE)/ni - ONE
  end do

  standard_op = StandardOperators1D(po_u)
  interpol_op = InterpolationOperator3D(standard_op, xi)

  allocate(vi(ni, ni, ni, ne, 3))
  allocate(wi(1:nl(1)  , 1:nl(2)  ,1:nl(3)  , 3))
  allocate(wt(0:nl(1)-1, 0:nl(2)-1,0:nl(3)/2, 3))

  plan = fftw_plan_dft_r2c_3d( ng(3), ng(2), ng(1), wi(:,:,:,1), wt(:,:,:,1) &
                             , FFTW_ESTIMATE )

  !!$omp parallel

  !-----------------------------------------------------------------------------
  ! FFT test

  if (any(kc >= 0) .or. any(ks > 0)) then

    ! global offsets
    i0 = nl(1) * ip(1) - 1
    j0 = nl(2) * ip(2) - 1
    k0 = nl(3) * ip(3) - 1

    do k = 1, nl(3)
    do j = 1, nl(2)
    do i = 1, nl(1)

      ig = i0 + i
      jg = j0 + j
      kg = k0 + k

      wi(i,j,k,1) = 0

      if (kc(1) >= 0)  wi(i,j,k,1) = wi(i,j,k,1) + cos(2*PI*kc(1)*ig/ng(1))
      if (kc(2) >= 0)  wi(i,j,k,1) = wi(i,j,k,1) + cos(2*PI*kc(2)*jg/ng(2))
      if (kc(3) >= 0)  wi(i,j,k,1) = wi(i,j,k,1) + cos(2*PI*kc(3)*kg/ng(3))

      if (ks(1) >  0)  wi(i,j,k,1) = wi(i,j,k,1) + sin(2*PI*ks(1)*ig/ng(1))
      if (ks(2) >  0)  wi(i,j,k,1) = wi(i,j,k,1) + sin(2*PI*ks(2)*jg/ng(2))
      if (ks(3) >  0)  wi(i,j,k,1) = wi(i,j,k,1) + sin(2*PI*ks(3)*kg/ng(3))

    end do
    end do
    end do

    call fftw_execute_dft_r2c(plan, wi(:,:,:,1), wt(:,:,:,1))

    do k = 0, ubound(wt,3)
    do j = 0, ubound(wt,2)
    do i = 0, ubound(wt,1)
      z = wt(i,j,k,1)
      if (abs(z) < epsilon(1.0)) cycle
      if (abs(z%re) < epsilon(1.0)) z%re = 0
      if (abs(z%im) < epsilon(1.0)) z%im = 0
      write(*,'(3I4,3X,2(2(ES12.5,2X),2X))') i, j, k, z, z / product(ng)
    end do
    end do
    end do

  end if

  ! interpolate to equidistant nodes
  call interpol_op % Apply(uc(:,:,:,:,1:3), vi)

  ! map to cartesian mesh
  !!$omp do collapse(3) private(e, i1, i2, j1, j2, k1, k2)
  do e3 = 1, ep(3)
  do e2 = 1, ep(2)
  do e1 = 1, ep(1)
    e  = e1 + ep(1) * (e2 - 1 + ep(2) * (e3 - 1))
    i1 = 1  + ep(1) * (e1 - 1)
    i2 = i1 + ni - 1
    j1 = 1  + ep(2) * (e2 - 1)
    j2 = j1 + ni - 1
    k1 = 1  + ep(3) * (e3 - 1)
    k2 = k1 + ni - 1
    wi(i1:i2,j1:j2,k1:k2,1:3) = vi(:,:,:,e,1:3)
  end do
  end do
  end do

  call fftw_execute_dft_r2c(plan, wi(:,:,:,1), wt(:,:,:,1))

  !!$omp end parallel

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!===============================================================================
! Standalone procedures

!-------------------------------------------------------------------------------
!> Loads the flow variables to disk

subroutine LoadFlowVariables(flow_case, t, u)
  character(len=*), intent(in)  :: flow_case    !< flow case
  real(RNP),        intent(out) :: t            !< time
  real(RNP),        intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: unit, ios
  character(len=100) :: file

  !$omp single

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'

  open(newunit=unit, file=file, status='OLD', form='UNFORMATTED', iostat=ios)
  if (ios == 0) then
    read(unit) t
    read(unit) u
    close(unit)
  else
    call Error('LoadFlowVariables', 'found no matching input file')
  end if

  !$omp end single

end subroutine LoadFlowVariables

!===============================================================================

end program Transition_TG_Post
