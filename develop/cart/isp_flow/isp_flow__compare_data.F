!> summary:  Comparison of incompressible single-phase flow data
!> author:   Joerg Stiller
!> date:     2021/05/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ISP_Flow__Compare_Data

  use Kind_Parameters
  use Execution_Control
  use Array_Reductions
  use XMPI

  implicit none

  ! declarations ...............................................................

  character(len=80) :: control_file   ! file containing control parameters
  character(len=80) :: ref_data_file  ! restart file with reference data
  character(len=80) :: tst_data_file  ! restart file with test data
  integer           :: po             ! polynomial order
  integer           :: n_elem         ! number of elements per partition

  ! parameters read from control file
  namelist /control/ ref_data_file, tst_data_file, po, n_elem

  real(RNP), allocatable :: u_ref(:,:,:,:,:)
  real(RNP), allocatable :: u_tst(:,:,:,:,:)

  type(MPI_Comm) :: comm
  integer        :: rank, n_proc, prm, stat
  real(RNP)      :: delta_v_loc, delta_v_max, delta_v_rms

  ! initialization .............................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  ! control parameters .........................................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Comparison of Cartesian IP/DG Incompressible Flow Data'
    write(*,*)

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'isp_flow__compare_data'
    end if

    write(*,'(2X,A,/)') 'reading ' // trim(control_file)//'.prm'
    open(newunit=prm, file=trim(control_file)//'.prm')
    read(prm, nml=control)
    close(prm)

    write(*,'(2X,A,T17,A,A)')  'ref_data_file', '=  ', trim(ref_data_file)
    write(*,'(2X,A,T17,A,A)')  'tst_data_file', '=  ', trim(tst_data_file)
    write(*,'(2X,A,T17,A,I0)') 'po           ', '=  ', po
    write(*,'(2X,A,T17,A,I0)') 'n_elem       ', '=  ', n_elem
    write(*,'(2X,A,T17,A,I0)') 'n_proc       ', '=  ', n_proc
    write(*,*)

  end if

  call XMPI_Bcast(ref_data_file, 0, comm)
  call XMPI_Bcast(tst_data_file, 0, comm)
  call XMPI_Bcast(po           , 0, comm)
  call XMPI_Bcast(n_elem       , 0, comm)

  ! load data ..................................................................

  allocate(u_ref(0:po,0:po,0:po,n_elem,4))
  allocate(u_tst(0:po,0:po,0:po,n_elem,4))
  call LoadFlowVariables(ref_data_file, u_ref)
  call LoadFlowVariables(tst_data_file, u_tst)

  ! compare data ...............................................................

  u_tst = u_tst - u_ref

  delta_v_loc = maxval(abs(u_tst))
  call XMPI_Reduce(delta_v_loc, delta_v_max, MPI_MAX, 0, comm)

  delta_v_rms = ScalarProduct(u_tst(:,:,:,:,1:3), u_tst(:,:,:,:,1:3), comm)
  delta_v_rms = sqrt( delta_v_rms / ((po+1)**3 * n_elem * n_proc) )

  ! print results ..............................................................

  if (rank == 0) then
    write(*,'(2X,A,T17,A,ES12.5)') 'delta_v_max', '= ', delta_v_max
    write(*,'(2X,A,T17,A,ES12.5)') 'delta_v_rms', '= ', delta_v_rms
    write(*,*)
  end if

  ! finalization ...............................................................

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Loads the flow variables

  subroutine LoadFlowVariables(data_file, u)
    character(len=*), intent(in)  :: data_file    !< restart file
    real(RNP),        intent(out) :: u(:,:,:,:,:) !< flow variables

    character(len=200) :: file
    real(RNP) :: t
    integer :: unit, ios

    write(file,'(2A,I0,A)') trim(data_file), '_p', rank, '.fld'
    open(newunit=unit, file=file, status='OLD', form='UNFORMATTED', iostat=ios)
    if (ios == 0) then
      read(unit) t
      read(unit) u
      close(unit)
    else
      call Error('SaveFlowVariables', 'found no matching input file')
    end if

  end subroutine LoadFlowVariables

  !=============================================================================

end program ISP_Flow__Compare_Data
