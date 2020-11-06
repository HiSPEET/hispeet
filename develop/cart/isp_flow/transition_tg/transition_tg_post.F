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

  type(MPI_Status) :: status

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

  ! global/transformed velocity on regular Cartesian mesh (rank 0 only)
  real    (C_DOUBLE)        , allocatable, save :: vg(:,:,:,:)
  complex (C_DOUBLE_COMPLEX), allocatable, save :: vt(:,:,:,:)

  logical   :: exists
  integer   :: prm, stat
  integer   :: e, e1, e2, e3
  integer   :: i, i0, ig
  integer   :: j, j0, jg
  integer   :: k, k0, kg

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
  call XMPI_Bcast(ks  , 0, comm)
  call XMPI_Bcast(kc  , 0, comm)

  ! derived parameters
  dx = (problem % x1 - problem % x0) / (np * ep)
  ne = product(ep)
  ni = po_u + 1
  nl = ni * ep
  ng = nl * np

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

  ! flow problem ...............................................................

  allocate(FlowProblem_Transition_TG :: problem)
  call problem % SetProblem(flow_case, comm)

  ! variables and operators ....................................................

  allocate(xi(ni))
  do i = 1, ni
    xi(i) = (TWO * i - ONE)/ni - ONE
  end do

  standard_op = StandardOperators1D(po_u)
  interpol_op = InterpolationOperator3D(standard_op, xi)

  allocate(vi(ni, ni, ni, ne, 3))

  if (rank == 0) then

    allocate(vg(1:ng(1)  , 1:ng(2)  , 1:ng(3)  , 3))
    allocate(vt(0:ng(1)/2, 0:ng(2)-1, 0:ng(3)-1, 3))

    plan = fftw_plan_dft_r2c_3d( ng(3), ng(2), ng(1), vg(:,:,:,1), vt(:,:,:,1) &
                               , FFTW_ESTIMATE )
  end if

  !-----------------------------------------------------------------------------
  ! Generate/read and interpolate flow data

  if (any(kc >= 0) .or. any(ks > 0)) then

    ! local FFT test data ......................................................

    ! global offsets
    i0 = nl(1) * ip(1)
    j0 = nl(2) * ip(2)
    k0 = nl(3) * ip(3)

    do e3 = 1, ep(3)
    do e2 = 1, ep(2)
    do e1 = 1, ep(1)

      do k = 1, ni
      do j = 1, ni
      do i = 1, ni

        ! element ID
        e  = e1 + ep(1) * (e2 - 1 + ep(2) * (e3 - 1))

        ! global point indices starting from 0
        ig = i0 + ni * (e1 - 1) + i - 1
        jg = j0 + ni * (e2 - 1) + j - 1
        kg = k0 + ni * (e3 - 1) + k - 1

        vi(i,j,k,e,1:3) = 0

        if (kc(1) >= 0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + cos(2*PI*kc(1)*ig/ng(1))
        if (kc(2) >= 0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + cos(2*PI*kc(2)*jg/ng(2))
        if (kc(3) >= 0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + cos(2*PI*kc(3)*kg/ng(3))

        if (ks(1) >  0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + sin(2*PI*ks(1)*ig/ng(1))
        if (ks(2) >  0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + sin(2*PI*ks(2)*jg/ng(2))
        if (ks(3) >  0) vi(i,j,k,e,1) = vi(i,j,k,e,1) + sin(2*PI*ks(3)*kg/ng(3))

      end do
      end do
      end do

    end do
    end do
    end do

  else

    ! local flow data ..........................................................

    ! load
    allocate(uc(0:po_u, 0:po_u, 0:po_u, ne, problem%nc))
    call LoadFlowVariables(flow_case, t, uc)

    ! interpolate to equidistant nodes
    call interpol_op % Apply(uc(:,:,:,:,1:3), vi)
    deallocate(uc)

  end if

  !-----------------------------------------------------------------------------
  ! Transfer to data to rank 0

  if (rank > 0) then

    call MPI_Send(vi, size(vi), MPI_REAL_RNP, 0, 0, MPI_COMM_WORLD)

  else
!print *, 'minval(vi(*,1)) = ', minval(vi(:,:,:,:,1))
!print *, 'maxval(vi(*,1)) = ', maxval(vi(:,:,:,:,1))

    ! own contribution
    call CopyLocalToGlobal(nl, ep, ip, vi, vg)

    do i = 1, product(np)-1

      ! receive data from partition i
      call MPI_Recv(vi, size(vi), MPI_REAL_RNP, i, 0, MPI_COMM_WORLD, status)

      ! copy into global velocity array
      ip = CartPID(i, np)
      call CopyLocalToGlobal(nl, ep, ip, vi, vg)

    end do

  end if

  !-----------------------------------------------------------------------------

!print *, 'minval(vg(*,1)) = ', minval(vg(:,:,:,1))
!print *, 'maxval(vg(*,1)) = ', maxval(vg(:,:,:,1))
  do i = 1, 3
    call fftw_execute_dft_r2c(plan, vg(:,:,:,i), vt(:,:,:,i))
  end do

  !-----------------------------------------------------------------------------
  ! Evaluation

  if (any(kc >= 0) .or. any(ks > 0)) then

    if (rank == 0) then
      write(*,'(/,A,/)') 'FFT test: active modes'
      write(*,'(  A  )') '   i   j   k       real          imag'
      do k = 0, ubound(vt,3)
      do j = 0, ubound(vt,2)
      do i = 0, ubound(vt,1)
        z = vt(i,j,k,1)
        if (abs(z) < epsilon(1.0)) cycle
        if (abs(z%re) < epsilon(1.0)) z%re = 0
        if (abs(z%im) < epsilon(1.0)) z%im = 0
        write(*,'(3I4,3X,2(ES12.5,2X))') i, j, k, z / product(ng)
      end do
      end do
      end do
    end if

  else

    ! TBD

  end if

  !-----------------------------------------------------------------------------
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
  !> Copy local velocity data to global cartesian array

  subroutine CopyLocalToGlobal(nl, ep, ip, vl, vg)
    integer,   intent(in)    :: nl(3)
    integer,   intent(in)    :: ep(3)
    integer,   intent(in)    :: ip(3)
    real(RNP), intent(in)    :: vl(:,:,:,:,:)
    real(RNP), intent(inout) :: vg(:,:,:,:)

    integer :: ni
    integer :: e1, e2, e3
    integer :: i0, il, ig
    integer :: j0, jl, jg
    integer :: k0, kl, kg

    ni = size(vl,1)

    ! global offsets
    i0 = nl(1) * ip(1)
    j0 = nl(2) * ip(2)
    k0 = nl(3) * ip(3)

    do e3 = 1, ep(3)
    do e2 = 1, ep(2)
    do e1 = 1, ep(1)

      do kl = 1, ni
      do jl = 1, ni
      do il = 1, ni

        ! element ID
        e  = e1 + ep(1) * (e2 - 1 + ep(2) * (e3 - 1))

        ! global point indices
        ig = i0 + ni * (e1 - 1) + il
        jg = j0 + ni * (e2 - 1) + jl
        kg = k0 + ni * (e3 - 1) + kl

        vg(ig,jg,kg,1:3) = vl(il,jl,kl,e,1:3)

      end do
      end do
      end do

    end do
    end do
    end do

  end subroutine CopyLocalToGlobal

  !-----------------------------------------------------------------------------
  !> Loads the flow variables from file

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

  !=============================================================================

end program Transition_TG_Post
