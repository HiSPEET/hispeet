module Data_Exchange__3D
  use Kind_Parameters
  use XMPI
  use Mesh_Map_To_Child__3D
  use Mesh_Map_To_Parent__3D

  implicit none
  private

  public :: DataExchangeMap_3D
  public :: DataExchangeSendBuf_3D
  public :: DataExchangeRecvBuf_3D

  !-----------------------------------------------------------------------------
  !> Map for exchanging element data

  type DataExchangeMap_3D
    type(MPI_Comm)       :: comm       !< MPI communicator
    integer              :: proc       !< MPI destination or source rank
    integer, allocatable :: id_elem(:) !< elements with data to be exchanged
  end type DataExchangeMap_3D

  ! constructors
  interface DataExchangeMap_3D
    module procedure New_ExchangeMap_from_MapToChild
    module procedure New_ExchangeMap_from_MapToParent
  end interface

  !-----------------------------------------------------------------------------
  !> Structure for sending data of retained elements to new location

  type DataExchangeSendBuf_3D
    class(DataExchangeMap_3D), pointer :: map => null()
    integer :: ne = -1 !< number of elements with data to be send
    integer :: na = -1 !< number of attribute entries per element
    integer :: np = -1 !< number of element points per direction
    integer :: nv = -1 !< number of variables
    integer,   allocatable :: buf_a(:,:)       !< send buffer for attributes
    real(RNP), allocatable :: buf_v(:,:,:,:,:) !< send buffer for variables
    type(MPI_Request)      :: req_a            !< MPI request for attributes
    type(MPI_Request)      :: req_v            !< MPI request for variables
  contains
    generic :: Extract_Data => Extract_ScalarData, Extract_ArrayData
    procedure, private :: Extract_ScalarData, Extract_ArrayData
    procedure :: Send_Start
    procedure :: Send_Finish
  end type DataExchangeSendBuf_3D

  !-----------------------------------------------------------------------------
  !> Structure for receiving data of retained elements from old location

  type DataExchangeRecvBuf_3D
    class(DataExchangeMap_3D), pointer :: map => null()
    integer :: ne = -1 !< number of elements with data to be received
    integer :: na = -1 !< number of attribute entries per element
    integer :: np = -1 !< number of element points per direction
    integer :: nv = -1 !< number of variables
    integer,   allocatable :: buf_a(:,:)       !< recv buffer for attributes
    real(RNP), allocatable :: buf_v(:,:,:,:,:) !< recv buffer for variables
    type(MPI_Request)      :: req_a            !< MPI request for attributes
    type(MPI_Request)      :: req_v            !< MPI request for variables
  contains
    generic :: Init => Init_Recv_ScalarData, Init_Recv_ArrayData
    procedure, private :: Init_Recv_ScalarData, Init_Recv_ArrayData
    procedure :: Recv_Start
    procedure :: Recv_Finish
    generic :: Assign_Data => Assign_ScalarData, Assign_ArrayData
    procedure, private :: Assign_ScalarData, Assign_ArrayData
  end type DataExchangeRecvBuf_3D

contains

  !=============================================================================
  ! TBP of DataExchangeMap_3D

  !-----------------------------------------------------------------------------
  !> New exchange map from map to child

  elemental function New_ExchangeMap_from_MapToChild(map_child, frozen) &
      result(this)

    class(MeshMapToChild_3D), intent(in) :: map_child
    logical, optional, intent(in) :: frozen !< in/exclude frozen elements [T]

    type(DataExchangeMap_3D) :: this
    integer :: n_elem

    ! default: send to all children
    n_elem = map_child % n_elem

    ! optional: exclude frozen children
    if (present(frozen)) then
      if (frozen) then
        n_elem = map_child % n_active
      end if
    end if

    this % comm    = map_child % comm
    this % proc    = map_child % proc
    this % id_elem = map_child % id_elem(:n_elem)

  end function New_ExchangeMap_from_MapToChild

  !-----------------------------------------------------------------------------
  !> New send map from map to parent

  elemental function New_ExchangeMap_from_MapToParent(map_parent, frozen) &
      result(this)

    class(MeshMapToParent_3D), intent(in) :: map_parent
    logical, optional, intent(in) :: frozen !< in/exclude frozen clusters [T]

    type(DataExchangeMap_3D) :: this
    integer :: n_cluster

    ! default: send to all parents
    n_cluster = map_parent % n_cluster

    ! optional: exclude frozen clusters
    if (present(frozen)) then
      if (frozen) then
        n_cluster = map_parent % n_active
      end if
    end if

    this % comm    = map_parent % comm
    this % proc    = map_parent % proc
    this % id_elem = map_parent % id_cluster(:n_cluster)

  end function New_ExchangeMap_from_MapToParent

  !=============================================================================
  ! TBP of DataExchangeSendBuf_3D

  !-----------------------------------------------------------------------------
  !> Extraction of attributes and/or scalar data for sending

  subroutine Extract_ScalarData(this, map, a, v)
    class(DataExchangeSendBuf_3D),     intent(inout) :: this
    class(DataExchangeMap_3D), target, intent(in)    :: map
    integer,     contiguous, optional, intent(in)    :: a(:,:)
    real(RNP),   contiguous, optional, intent(in)    :: v(:,:,:,:)

    this % map => map
    this % ne  =  size(map % id_elem)

    call Extract_Attributes(this, a)
    call Extract_ScalarVariable(this, v)

  end subroutine Extract_ScalarData

  !-----------------------------------------------------------------------------
  !> Extraction of attributes and/or array data for sending

  subroutine Extract_ArrayData(this, map, a, v)
    class(DataExchangeSendBuf_3D),     intent(inout) :: this
    class(DataExchangeMap_3D), target, intent(in)    :: map
    integer,     contiguous, optional, intent(in)    :: a(:,:)
    real(RNP),   contiguous,           intent(in)    :: v(:,:,:,:,:)

    this % map => map
    this % ne  =  size(map % id_elem)

    call Extract_Attributes(this, a)
    call Extract_ArrayVariable(this, v)

  end subroutine Extract_ArrayData

  !-----------------------------------------------------------------------------
  !> Extraction of attributes for sending

  subroutine Extract_Attributes(this, a)
    class(DataExchangeSendBuf_3D), intent(inout) :: this
    integer, contiguous, optional, intent(in) :: a(:,:)

    integer :: e

    associate( ne => this % ne &
             , na => this % na &
             , id_elem => this % map % id_elem )

      if (present(a)) then

        na = size(a,1)

        ! check and possibly (re)allocate attribute send buffer
        if (allocated(this % buf_a)) then
          if (any(shape(this % buf_a) /= [na,ne])) then
            deallocate(this % buf_a)
          end if
        end if
        if (.not. allocated(this % buf_a)) then
          allocate(this % buf_a(na,ne))
        end if

        ! extract attributes to send buffer
        do e = 1, ne
          this % buf_a(:,e) = a(:,id_elem(e))
        end do

      else

        na = 0
        if (allocated(this % buf_a)) then
          deallocate(this % buf_a)
        end if

      end if
    end associate

  end subroutine Extract_Attributes

  !-----------------------------------------------------------------------------
  !> Extraction of scalar data for sending

  subroutine Extract_ScalarVariable(this, v)
    class(DataExchangeSendBuf_3D), intent(inout) :: this
    real(RNP), contiguous, optional, intent(in) :: v(:,:,:,:)

    integer :: e

    associate( ne => this % ne &
             , np => this % np &
             , nv => this % nv &
             , id_elem => this % map % id_elem )

      if (present(v)) then

        np = size(v,1)
        nv = 1

        ! check and possibly (re)allocate variable send buffer
        if (allocated(this % buf_v)) then
          if (any(shape(this % buf_v) /= [np,np,np,ne,nv])) then
            deallocate(this % buf_v)
          end if
        end if
        if (.not. allocated(this % buf_v)) then
          allocate(this % buf_v(np,np,np,ne,nv))
        end if

        ! extract variables to send buffer
        do e = 1, ne
          this % buf_v(:,:,:,e,1) = v(:,:,:,id_elem(e))
        end do

      else

        np = 0
        nv = 0
        if (allocated(this % buf_v)) then
          deallocate(this % buf_v)
        end if

      end if
    end associate

  end subroutine Extract_ScalarVariable

  !-----------------------------------------------------------------------------
  !> Extraction of array data for sending

  subroutine Extract_ArrayVariable(this, v)
    class(DataExchangeSendBuf_3D), intent(inout) :: this
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)

    integer :: e, k

    associate( ne => this % ne &
             , np => this % np &
             , nv => this % nv &
             , id_elem => this % map % id_elem )

      np = size(v,1)
      nv = size(v,5)

      ! check and possibly (re)allocate variable send buffer
      if (allocated(this % buf_v)) then
        if (any(shape(this % buf_v) /= [np,np,np,ne,nv])) then
          deallocate(this % buf_v)
        end if
      end if
      if (.not. allocated(this % buf_v)) then
        allocate(this % buf_v(np,np,np,ne,nv))
      end if

      ! extract variables to send buffer
      do k = 1, nv
      do e = 1, ne
        this % buf_v(:,:,:,e,k) = v(:,:,:,id_elem(e),k)
      end do
      end do

    end associate

  end subroutine Extract_ArrayVariable

  !-----------------------------------------------------------------------------
  !> Nonblocking sending of attributes and data

  subroutine Send_Start(this)
    class(DataExchangeSendBuf_3D), asynchronous, intent(inout) :: this

    if (this % na > 0) then
      call XMPI_Isend(this%buf_a, this%map%proc, 71, this%map%comm, this%req_a)
    end if

    if (this % nv > 0) then
      call XMPI_Isend(this%buf_v, this%map%proc, 72, this%map%comm, this%req_v)
    end if

  end subroutine Send_Start

  !-----------------------------------------------------------------------------
  !> Finalization of sending

  subroutine Send_Finish(this)
    class(DataExchangeSendBuf_3D), asynchronous, intent(inout) :: this

    if (this % na > 0) call MPI_Wait(this % req_a, MPI_STATUS_IGNORE)
    if (this % nv > 0) call MPI_Wait(this % req_v, MPI_STATUS_IGNORE)

  end subroutine Send_Finish

  !=============================================================================
  ! TBP of DataExchangeRecvBuf_3D

  !-----------------------------------------------------------------------------
  !> Initialization of receive buffer for attributes and/or scalar data

  subroutine Init_Recv_ScalarData(this, map, a, v)
    class(DataExchangeRecvBuf_3D),     intent(inout) :: this
    class(DataExchangeMap_3D), target, intent(in)    :: map
    integer,     contiguous, optional, intent(in)    :: a(:,:)
    real(RNP),   contiguous, optional, intent(in)    :: v(:,:,:,:)

    this % map => map
    this % ne  =  size(map % id_elem)

    call Init_Recv_Attributes(this, a)
    call Init_Recv_ScalarVariable(this, v)

  end subroutine Init_Recv_ScalarData

  !-----------------------------------------------------------------------------
  !> Initialization of receive buffer for attributes and/or array data

  subroutine Init_Recv_ArrayData(this, map, a, v)
    class(DataExchangeRecvBuf_3D),         intent(inout) :: this
    class(DataExchangeMap_3D), target, intent(in)    :: map
    integer,         contiguous, optional, intent(in)    :: a(:,:)
    real(RNP),       contiguous,           intent(in)    :: v(:,:,:,:,:)

    this % map => map
    this % ne  =  size(map % id_elem)

    call Init_Recv_Attributes(this, a)
    call Init_Recv_ArrayVariable(this, v)

  end subroutine Init_Recv_ArrayData

  !-----------------------------------------------------------------------------
  !> Initialization of receive buffer for attributes

  subroutine Init_Recv_Attributes(this, a)
    class(DataExchangeRecvBuf_3D), intent(inout) :: this
    integer, contiguous, optional, intent(in) :: a(:,:)

    associate( ne => this % ne &
             , na => this % na )

      if (present(a)) then

        na = size(a,1)

        if (allocated(this % buf_a)) then
          if (any(shape(this % buf_a) /= [na,ne])) then
            deallocate(this % buf_a)
          end if
        end if
        if (.not. allocated(this % buf_a)) then
          allocate(this % buf_a(na,ne))
        end if

      else

        na = 0
        if (allocated(this % buf_a)) then
          deallocate(this % buf_a)
        end if

      end if
    end associate

  end subroutine Init_Recv_Attributes

  !-----------------------------------------------------------------------------
  !> Initialization of receive buffer for scalar data

  subroutine Init_Recv_ScalarVariable(this, v)
    class(DataExchangeRecvBuf_3D), intent(inout) :: this
    real(RNP), contiguous, optional, intent(in) :: v(:,:,:,:)

    associate( ne => this % ne &
             , np => this % np &
             , nv => this % nv )

      if (present(v)) then

        np = size(v,1)
        nv = 1

        if (allocated(this % buf_v)) then
          if (any(shape(this % buf_v) /= [np,np,np,ne,nv])) then
            deallocate(this % buf_v)
          end if
        end if
        if (.not. allocated(this % buf_v)) then
          allocate(this % buf_v(np,np,np,ne,nv))
        end if

      else

        np = 0
        nv = 0
        if (allocated(this % buf_v)) then
          deallocate(this % buf_v)
        end if

      end if
    end associate

  end subroutine Init_Recv_ScalarVariable

  !-----------------------------------------------------------------------------
  !> Initialization of receive buffer for array data

  subroutine Init_Recv_ArrayVariable(this, v)
    class(DataExchangeRecvBuf_3D), intent(inout) :: this
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)

    associate( ne => this % ne &
             , np => this % np &
             , nv => this % nv )

      np = size(v,1)
      nv = size(v,5)

      if (allocated(this % buf_v)) then
        if (any(shape(this % buf_v) /= [np,np,np,ne,nv])) then
          deallocate(this % buf_v)
        end if
      end if
      if (.not. allocated(this % buf_v)) then
        allocate(this % buf_v(np,np,np,ne,nv))
      end if

    end associate

  end subroutine Init_Recv_ArrayVariable

  !-----------------------------------------------------------------------------
  !> Nonblocking receiving of attributes and data

  subroutine Recv_Start(this)
    class(DataExchangeRecvBuf_3D), asynchronous, intent(inout) :: this

    if (this % na > 0) then
      call XMPI_Irecv(this%buf_a, this%map%proc, 71, this%map%comm, this%req_a)
    end if

    if (this % nv > 0) then
      call XMPI_Irecv(this%buf_v, this%map%proc, 72, this%map%comm, this%req_v)
    end if

  end subroutine Recv_Start

  !-----------------------------------------------------------------------------
  !> Finalization of receiving

  subroutine Recv_Finish(this)
    class(DataExchangeRecvBuf_3D), asynchronous, intent(inout) :: this

    if (this % na > 0) call MPI_Wait(this % req_a, MPI_STATUS_IGNORE)
    if (this % nv > 0) call MPI_Wait(this % req_v, MPI_STATUS_IGNORE)

  end subroutine Recv_Finish

  !-----------------------------------------------------------------------------
  !> Assignment of received attributes and/or scalar data

  subroutine Assign_ScalarData(this, a, v)
    class(DataExchangeRecvBuf_3D), asynchronous, intent(in) :: this
    integer,   contiguous, optional, intent(inout) :: a(:,:)
    real(RNP), contiguous, optional, intent(inout) :: v(:,:,:,:)

    integer :: e

    associate( ne => this % ne &
             , na => this % na &
             , np => this % np &
             , nv => this % nv &
             , id_elem => this % map % id_elem )

      if (present(a) .and. na > 0) then
        do e = 1, ne
          a(:,id_elem(e)) = this % buf_a(:,e)
        end do
      end if

      if (present(v) .and. nv > 0) then
        do e = 1, ne
          v(:,:,:,id_elem(e)) = this % buf_v(:,:,:,e,1)
        end do
      end if

    end associate

  end subroutine Assign_ScalarData

  !-----------------------------------------------------------------------------
  !> Assignment of received attributes and/or array data

  subroutine Assign_ArrayData(this, a, v)
    class(DataExchangeRecvBuf_3D), asynchronous, intent(in) :: this
    integer,   contiguous, optional, intent(inout) :: a(:,:)
    real(RNP), contiguous,           intent(inout) :: v(:,:,:,:,:)

    integer :: e, k

    associate( ne => this % ne &
             , na => this % na &
             , np => this % np &
             , nv => this % nv &
             , id_elem => this % map % id_elem )

      if (present(a) .and. na > 0) then
        do e = 1, ne
          a(:,id_elem(e)) = this % buf_a(:,e)
        end do
      end if

      do k = 1, nv
      do e = 1, ne
        v(:,:,:,id_elem(e),k) = this % buf_v(:,:,:,e,k)
      end do
      end do

    end associate

  end subroutine Assign_ArrayData

  !=============================================================================

end module Data_Exchange__3D
