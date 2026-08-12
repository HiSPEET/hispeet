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

!> summary:  Change of polynomial order
!>
!> Loads restart files and changes the polynom order. If present, the first
!> argument of the invoking command will be interpreted as the base name of the
!> control file. If absent, the program looks for `change_polynomial_order.prm`.
!===============================================================================

program Change_Polynomial_Order

  use, intrinsic :: IEEE_Arithmetic, only: ieee_is_nan

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use OpenMP_Binding
  use Execution_Control
  use Array_Assignments
  use Array_Reductions
  use Gauss_Jacobi
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__3D
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
  character(len=80) :: flow_case    = ''
  character(len=80) :: solver_file  = 'isp_flow__sdc_test'
  logical           :: reset_time   = .false. ! resets time to zero
  logical           :: delete_avg   = .false. ! deletes time-average of input data
  integer           :: po_u_new     =  5      ! new polinomal order

  ! parameters read from control file
  namelist /control/ flow_case, solver_file, reset_time, delete_avg, po_u_new

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
  integer   :: ne           ! number of elements in the local partition

  ! flow problem ...............................................................

  class(FlowProblem), allocatable :: problem

  ! operators ..................................................................

  type(StandardElementOperators_1D)    , save :: standard_op
  type(StandardElementOperators_1D)    , save :: standard_op_new
  type(EmbeddedInterpolationOperator_3D), save :: interpol_op

  ! variables ..................................................................

  ! computed solution from input flow data
  real(RNP), allocatable, save :: u(:,:,:,:,:)        ! computed solution
  real(RNP), allocatable, save :: q_avg(:,:,:,:,:)    ! time-averaged quantities
  real(RNP)                    :: t
  integer                      :: n_avg

  ! interpolated solution
  real(RNP), allocatable, save :: ui(:,:,:,:,:)
  real(RNP), allocatable, save :: q_avg_i(:,:,:,:,:)  ! interpolated time-averaged quantities
  real(RNP), allocatable, save :: xi(:)               ! interpolation points in [-1,1]

  logical   :: exists
  integer   :: prm, stat
  integer   :: i

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
    write(*,'(A)') 'Changing Polynomial Order'
    write(*,*)

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'change_polynomial_order'
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
    write(*,'(2X,A,1X,I0)') 'po_old      = ', po_u
    write(*,'(2X,A,1X,I0)') 'po_new      = ', po_u_new
    write(*,*)

  end if

  ! globalize parameters
  call XMPI_Bcast(flow_case , 0, comm)
  call XMPI_Bcast(reset_time, 0, comm)
  call XMPI_Bcast(delete_avg, 0, comm)
  call XMPI_Bcast(ep        , 0, comm)
  call XMPI_Bcast(po_u      , 0, comm)
  call XMPI_Bcast(po_u_new  , 0, comm)

  ! flow problem ...............................................................

  allocate(FlowProblem_Channel :: problem)

  call problem % SetProblem(flow_case, comm)

  ! derived parameters .........................................................

  ne = product(ep)

  ! variables and operators ....................................................

  standard_op_new = StandardElementOperators_1D(po_u_new)

  allocate(xi(po_u_new+1))
  do i = 1, po_u_new+1
    xi(i) = standard_op_new%x(i-1)
  end do

  standard_op = StandardElementOperators_1D(po_u)
  interpol_op = EmbeddedInterpolationOperator_3D(standard_op, xi)

  !-----------------------------------------------------------------------------
  ! Read, interpolate and save flow data

  ! load
  if (rank == 0) then
    write(*,'(A)') 'Loading restart data'
  end if

  allocate(     u(0:po_u, 0:po_u, 0:po_u, ne, problem%nc),  &
            q_avg(0:po_u, 0:po_u, 0:po_u, ne, 9)            )

  call LoadFlowVariables(t, u, n_avg, q_avg)

  ! interpolate
  if (rank == 0) then
    write(*,'(A)') 'Interpolation to new polynomal order'
  end if

  allocate( ui(0:po_u_new, 0:po_u_new, 0:po_u_new, ne, problem%nc))
  call interpol_op % Apply(u(:,:,:,:,:), ui)

  if (n_avg > 0) then
    if (delete_avg) then
      q_avg = 0
      n_avg = 0
    else
    call interpol_op % Apply(q_avg(:,:,:,:,:), q_avg_i)
    end if
  end if

  if (reset_time) then
    t = 0
  end if

  ! save
  if (rank == 0) then
    write(*,'(A)') 'Save restart data'
  end if

  call SaveFlowVariables(t, u, n_avg, q_avg)

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> Loads the flow variables to disk

subroutine LoadFlowVariables(t, u, n_avg, q_avg)
  real(RNP),           intent(out) :: t                    !< time
  real(RNP),           intent(out) :: u(:,:,:,:,:)         !< flow variables
  integer  ,           intent(out) :: n_avg                !< current num averages
  real(RNP), optional, intent(out) :: q_avg(:,:,:,:,:)     !< time-averaged quantities

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
      if (n_avg > 0) then
        read(unit) q_avg
      else
        call SetArray(q_avg, ZERO, multi = .true.)
      end if
    end if

    close(unit)
  else
    call Error('SaveFlowVariables', 'found no matching input file')
  end if

  !$omp end single

end subroutine LoadFlowVariables

!-------------------------------------------------------------------------------
!> Saves the flow variables to disk

subroutine SaveFlowVariables(t, u, n_avg, q_avg)
  real(RNP),           intent(in) :: t                     !< time
  real(RNP),           intent(in) :: u(:,:,:,:,:)          !< flow variables
  integer  ,           intent(in) :: n_avg                 !< current num averages
  real(RNP), optional, intent(in) :: q_avg(:,:,:,:,:)      !< time-averaged quantities

  integer :: unit
  character(len=100) :: file

  !$omp single

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'
  open(newunit=unit, file=file, status='REPLACE', form='UNFORMATTED')
  write(unit) t
  write(unit) u
  if (present(q_avg)) then
    write(unit) n_avg
    write(unit) q_avg
  else
    write(unit) 0
  end if

  close(unit)

  !$omp end single

end subroutine SaveFlowVariables

!=============================================================================

end program Change_Polynomial_Order
