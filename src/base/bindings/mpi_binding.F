!> summary:  Initialization of MPI and provison of custom intrinsic types
!> author:   Joerg Stiller
!> date:     2014/01/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   * The real MPI type should be generated with `MPI_Type_create_f90_real`.
!>     As this fails with Intel's MPI, a workaround must be used :(
!==============================================================================

module MPI_Binding
  use MPI_F08
  use Kind_Parameters, only: IXS, IXL, RSP, RDP, RHP, RNP
  public

  private :: IXS, IXL, RSP, RDP, RNP, RHP

  !-----------------------------------------------------------------------------
  ! MPI data types

  type(MPI_Datatype), save :: &
    MPI_INTEGER_IXS, & !< MPI data type matching integer(IXS)
    MPI_INTEGER_IXL, & !< MPI data type matching integer(IXL)
    MPI_REAL_RSP,    & !< MPI data type matching real(RDP)
    MPI_REAL_RDP,    & !< MPI data type matching real(RDP)
    MPI_REAL_RHP,    & !< MPI data type matching real(RHP)
    MPI_REAL_RNP       !< MPI data type matching real(RNP)

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

  !-----------------------------------------------------------------------------
  !> Initializes the MPI data types defined in the module
  !>
  !> Initializes the MPI environment via MPI_Init_Thread and creates the
  !> datatypes used for MPI communication. If OpenMP is utilized, the MPI
  !> implementation is required to provide threading level MPI_THREAD_FUNNELED
  !> or higher.

  subroutine Init_MPI_Binding()
    integer :: provided_threading ! threading support, MPI 3.0 p.485-488

    ! safeguard ................................................................

    if (initialized) return

    call MPI_Initialized(initialized)

    ! initialize MPI with required threading support ...........................

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

    ! integer types ............................................................

    call MPI_Type_create_f90_integer(range(1_IXS), MPI_INTEGER_IXS)
    call MPI_Type_create_f90_integer(range(1_IXL), MPI_INTEGER_IXL)

    ! real types ...............................................................

    MPI_REAL_RSP = MPI_REAL
    MPI_REAL_RDP = MPI_DOUBLE_PRECISION

    select case(RHP)
    case(RDP)
      MPI_REAL_RHP = MPI_REAL_RDP
    case default
      MPI_REAL_RHP = MPI_REAL16
    end select

    select case(RNP)
    case(RSP)
      MPI_REAL_RNP = MPI_REAL_RSP
    case(RDP)
      MPI_REAL_RNP = MPI_REAL_RDP
    case(RHP)
      MPI_REAL_RNP = MPI_REAL_RHP
    end select

    initialized = .true.

  end subroutine Init_MPI_Binding

  !============================================================================

end module MPI_Binding
