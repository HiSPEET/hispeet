!> summary:  Extended MPI Fortran binding for type character
!> author:   Joerg Stiller
!> date:     2015/05/09, revised 2017/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Extended MPI Fortran binding for type character
!===============================================================================

module XMPI__Character
  use MPI_Binding
  implicit none
  private

  public :: XMPI_Bcast
  public :: XMPI_Ibcast
  public :: XMPI_Isend
  public :: XMPI_Irecv

  !-----------------------------------------------------------------------------
  !> Extended MPI_Bcast for type character

  interface XMPI_Bcast
    module procedure BcastX0
    module procedure BcastX1
    module procedure BcastX2
  end interface XMPI_Bcast

  !-----------------------------------------------------------------------------
  !> Extended MPI_Ibcast for type character

  interface XMPI_Ibcast
    module procedure IbcastX0
    module procedure IbcastX1
    module procedure IbcastX2
  end interface XMPI_Ibcast

  !-----------------------------------------------------------------------------
  !> Enhanced MPI_Isend for type character

  interface XMPI_Isend
    module procedure IsendX0
    module procedure IsendX1
    module procedure IsendX2
    module procedure IsendX3
    module procedure IsendX4
  end interface XMPI_Isend

  !-----------------------------------------------------------------------------
  !> Enhanced MPI_Irecv for type character

  interface XMPI_Irecv
    module procedure IrecvX0
    module procedure IrecvX1
    module procedure IrecvX2
    module procedure IrecvX3
    module procedure IrecvX4
  end interface XMPI_Irecv

contains

!===============================================================================
! XMPI_Bcast

!-------------------------------------------------------------------------------
!> Bcast for scalar buffers

subroutine BcastX0(buffer, root, comm)
  character(len=*), intent(inout) :: buffer !< buffer
  integer,          intent(in)    :: root   !< rank of broadcast root
  type(MPI_Comm),   intent(in)    :: comm   !< communicator

  call MPI_Bcast(buffer, len(buffer), MPI_CHARACTER, root, comm)

end subroutine BcastX0

!-------------------------------------------------------------------------------
!> Bcast for 1D buffers

subroutine BcastX1(buffer, root, comm)
  character(len=*), intent(inout) :: buffer(:) !< buffer
  integer,          intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm),   intent(in)    :: comm      !< communicator

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Bcast(buffer, l, MPI_CHARACTER, root, comm)

end subroutine BcastX1

!-------------------------------------------------------------------------------
!> Bcast for 2D buffers

subroutine BcastX2(buffer, root, comm)
  character(len=*), intent(inout) :: buffer(:,:) !< buffer
  integer,          intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm),   intent(in)    :: comm        !< communicator

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Bcast(buffer, l, MPI_CHARACTER, root, comm)

end subroutine BcastX2

!===============================================================================
! XMPI_Ibcast

!-------------------------------------------------------------------------------
!> Ibcast for scalar buffers

subroutine IbcastX0(buffer, root, comm, request)
  character(len=*),  intent(inout) :: buffer   !< buffer
  integer,           intent(in)    :: root     !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm     !< communicator
  type(MPI_Request), intent(out)   :: request  !< request

  call MPI_Ibcast(buffer, len(buffer), MPI_CHARACTER, root, comm, request)

end subroutine IbcastX0

!-------------------------------------------------------------------------------
!> Ibcast for 1D buffers

subroutine IbcastX1(buffer, root, comm, request)
  character(len=*),  intent(inout) :: buffer(:) !< buffer
  integer,           intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm      !< communicator
  type(MPI_Request), intent(out)   :: request   !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Ibcast(buffer, l, MPI_CHARACTER, root, comm, request)

end subroutine IbcastX1

!-------------------------------------------------------------------------------
!> Ibcast for 2D buffers

subroutine IbcastX2(buffer, root, comm, request)
  character(len=*),  intent(inout) :: buffer(:,:) !< buffer
  integer,           intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm        !< communicator
  type(MPI_Request), intent(out)   :: request     !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Ibcast(buffer, l, MPI_CHARACTER, root, comm, request)

end subroutine IbcastX2

!===============================================================================
! XMPI_Isend

!-------------------------------------------------------------------------------
!> Isend for scalar buffers

subroutine IsendX0(buffer, dest, tag, comm, request)
  character(len=*),  intent(in)  :: buffer   !< send buffer
  integer,           intent(in)  :: dest     !< rank of the receiver
  integer,           intent(in)  :: tag      !< message tag
  type(MPI_Comm),    intent(in)  :: comm     !< communicator
  type(MPI_Request), intent(out) :: request  !< request

  integer :: l

  l = len(buffer)
  call MPI_Isend(buffer, l, MPI_CHARACTER, dest, tag, comm, request)

end subroutine IsendX0

!-------------------------------------------------------------------------------
!> Isend for 1D buffers

subroutine IsendX1(buffer, dest, tag, comm, request)
  character(len=*),  intent(in)  :: buffer(:) !< send buffer
  integer,           intent(in)  :: dest      !< rank of the receiver
  integer,           intent(in)  :: tag       !< message tag
  type(MPI_Comm),    intent(in)  :: comm      !< communicator
  type(MPI_Request), intent(out) :: request   !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Isend(buffer, l, MPI_CHARACTER, dest, tag, comm, request)

end subroutine IsendX1

!-------------------------------------------------------------------------------
!> Isend for 2D buffers

subroutine IsendX2(buffer, dest, tag, comm, request)
  character(len=*),  intent(in)  :: buffer(:,:) !< send buffer
  integer,           intent(in)  :: dest        !< rank of the receiver
  integer,           intent(in)  :: tag         !< message tag
  type(MPI_Comm),    intent(in)  :: comm        !< communicator
  type(MPI_Request), intent(out) :: request     !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Isend(buffer, l, MPI_CHARACTER, dest, tag, comm, request)

end subroutine IsendX2

!-------------------------------------------------------------------------------
!> MPI_Isend for 3D buffers

subroutine IsendX3(buffer, dest, tag, comm, request)
  character(len=*),  intent(in)  :: buffer(:,:,:) !< send buffer
  integer,           intent(in)  :: dest          !< rank of the receiver
  integer,           intent(in)  :: tag           !< message tag
  type(MPI_Comm),    intent(in)  :: comm          !< communicator
  type(MPI_Request), intent(out) :: request       !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Isend(buffer, l, MPI_CHARACTER, dest, tag, comm, request)

end subroutine IsendX3

!-------------------------------------------------------------------------------
!> Isend for 4D buffers

subroutine IsendX4(buffer, dest, tag, comm, request)
  character(len=*),  intent(in)  :: buffer(:,:,:,:) !< send buffer
  integer,           intent(in)  :: dest            !< rank of the receiver
  integer,           intent(in)  :: tag             !< message tag
  type(MPI_Comm),    intent(in)  :: comm            !< communicator
  type(MPI_Request), intent(out) :: request         !< request

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Isend(buffer, l, MPI_CHARACTER, dest, tag, comm, request)

end subroutine IsendX4

!===============================================================================
! XMPI_Irecv

!-------------------------------------------------------------------------------
!> Irecv for scalar buffers

subroutine IrecvX0(buffer, source, tag, comm, request)
  character(len=*),  intent(out) :: buffer   !< reveive buffer
  integer,           intent(in)  :: source   !< rank of the sender
  integer,           intent(in)  :: tag      !< message tag
  type(MPI_Comm),    intent(in)  :: comm     !< communicator
  type(MPI_Request), intent(out) :: request  !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  integer :: l

  l = len(buffer)
  call MPI_Irecv(buffer, l, MPI_CHARACTER, source, tag, comm, request)

end subroutine IrecvX0

!-------------------------------------------------------------------------------
!> Irecv for 1D buffers

subroutine IrecvX1(buffer, source, tag, comm, request)
  character(len=*),  intent(out) :: buffer(:) !< reveive buffer
  integer,           intent(in)  :: source    !< rank of the sender
  integer,           intent(in)  :: tag       !< message tag
  type(MPI_Comm),    intent(in)  :: comm      !< communicator
  type(MPI_Request), intent(out) :: request   !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Irecv(buffer, l, MPI_CHARACTER, source, tag, comm, request)

end subroutine IrecvX1

!-------------------------------------------------------------------------------
!> Irecv for 2D buffers

subroutine IrecvX2(buffer, source, tag, comm, request)
  character(len=*),  intent(out) :: buffer(:,:) !< reveive buffer
  integer,           intent(in)  :: source      !< rank of the sender
  integer,           intent(in)  :: tag         !< message tag
  type(MPI_Comm),    intent(in)  :: comm        !< communicator
  type(MPI_Request), intent(out) :: request     !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Irecv(buffer, l, MPI_CHARACTER, source, tag, comm, request)

end subroutine IrecvX2

!-------------------------------------------------------------------------------
!> Irecv for 3D buffers

subroutine IrecvX3(buffer, source, tag, comm, request)
  character(len=*),  intent(out) :: buffer(:,:,:) !< reveive buffer
  integer,           intent(in)  :: source        !< rank of the sender
  integer,           intent(in)  :: tag           !< message tag
  type(MPI_Comm),    intent(in)  :: comm          !< communicator
  type(MPI_Request), intent(out) :: request       !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Irecv(buffer, l, MPI_CHARACTER, source, tag, comm, request)

end subroutine IrecvX3

!-------------------------------------------------------------------------------
!> Irecv for 4D buffers

subroutine IrecvX4(buffer, source, tag, comm, request)
  character(len=*),  intent(out) :: buffer(:,:,:,:) !< reveive buffer
  integer,           intent(in)  :: source          !< rank of the sender
  integer,           intent(in)  :: tag             !< message tag
  type(MPI_Comm),    intent(in)  :: comm            !< communicator
  type(MPI_Request), intent(out) :: request         !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  integer :: l

  l = len(buffer) * size(buffer)
  call MPI_Irecv(buffer, l, MPI_CHARACTER, source, tag, comm, request)

end subroutine IrecvX4

!===============================================================================

end module XMPI__Character
