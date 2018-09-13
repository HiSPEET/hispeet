!> summary:  Initialization of MPI and provison of custom intrinsic types
!> author:   Joerg Stiller
!> date:     2014/01/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Initialization of MPI and provison of custom intrinsic types
!==============================================================================

module MPI_Binding
  use MPI_F08
  use Kind_Parameters, only: IXS, IXL, RNP, RHP, RDP
  public

  private :: IXS, IXL, RNP, RHP, RDP

  !-----------------------------------------------------------------------------
  ! MPI data types

  type(MPI_Datatype), save :: &
    MPI_INTEGER_IXS, & !< MPI data type matching integer(IXS)
    MPI_INTEGER_IXL, & !< MPI data type matching integer(IXL)
    MPI_REAL_RDP,    & !< MPI data type matching real(RDP)
    MPI_REAL_RNP,    & !< MPI data type matching real(RNP)
    MPI_REAL_RHP       !< MPI data type matching real(RHP)

  !-----------------------------------------------------------------------------
  ! Private variables

  logical, private :: initialized = .false.

#ifdef _OPENMP
  ! Only the master thread may execute MPI routines
  integer, parameter :: REQUIRED_THREADING = MPI_THREAD_FUNNELED
#else
  ! Only one thread present
  integer, parameter :: REQUIRED_THREADING = MPI_THREAD_SINGLE
#endif

contains

!-------------------------------------------------------------------------------
!> Initializes the MPI data types defined in the module
!>
!> Initializes the MPI environment via MPI_Init_Thread and creates the datatypes
!> used for MPI communication. If OpenMP is utilized, the MPI implementation is
!> required to provide threading level MPI_THREAD_FUNNELED or higher.

subroutine Init_MPI_Binding()
  integer :: provided_threading ! threading support, MPI 3.0 p.485-488

  ! safeguard ..................................................................

  if (initialized) return

  call MPI_Initialized(initialized)

  ! initialize MPI with required threading support .............................

  ! when using OpenMP *only* the master thread may execute MPI calls.
  if(.not. initialized) then
    call MPI_Init_Thread(REQUIRED_THREADING, provided_threading)
  else
    call MPI_Query_thread(provided_threading)
  end if

  ! bail out if requirements are not met
  ! MPI_THREAD_SINGLE < MPI_THREAD_FUNNELED < MPI_THREAD_SERIALIZED
  if (provided_threading < REQUIRED_THREADING) then
    call MPI_Finalize()
    stop 'Init_MPI_Binding: insufficient threading support, aborting'
  end if

  ! integer types ..............................................................

  call MPI_Type_create_f90_integer(range(1_IXS), MPI_INTEGER_IXS)
  call MPI_Type_create_f90_integer(range(1_IXL), MPI_INTEGER_IXL)

  ! real types .................................................................

  MPI_REAL_RDP = MPI_Real_Type(RDP, range(1.0_RDP))
  MPI_REAL_RNP = MPI_Real_Type(RNP, range(1.0_RNP))
  MPI_REAL_RHP = MPI_Real_Type(RHP, range(1.0_RHP))

  initialized = .true.

contains

  !> Automatic selection/creation of matching MPI real datatype
  type(MPI_Datatype) function MPI_Real_Type(k, r)
    integer, intent(in) :: k  !< kind parameter
    integer, intent(in) :: r  !< precision

    if (k == kind(1E0)) then
      MPI_Real_Type = MPI_REAL
    else if (k == kind(1D0)) then
      MPI_Real_Type = MPI_DOUBLE_PRECISION
    else
      call MPI_Type_create_f90_real(MPI_UNDEFINED, r, MPI_Real_Type)
    end if

  end function MPI_Real_Type

end subroutine Init_MPI_Binding

!==============================================================================

end module MPI_Binding
