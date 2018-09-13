!> summary:  Extended MPI Fortran binding for type real(RNP)
!> author:   Joerg Stiller
!> date:     2014/11/06, revised 2017/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Extended MPI Fortran binding for type real(RNP)
!===============================================================================

module XMPI__Real_RNP
  use Kind_Parameters, only: RNP
  use MPI_Binding
  implicit none
  private

  public :: XMPI_Bcast
  public :: XMPI_Ibcast
  public :: XMPI_Isend
  public :: XMPI_Irecv
  public :: XMPI_Reduce
  public :: XMPI_Allreduce

  !-----------------------------------------------------------------------------
  !> Extended MPI_Bcast for type real(RNP)

  interface XMPI_Bcast
    module procedure BcastX0
    module procedure BcastX1
    module procedure BcastX2
  end interface XMPI_Bcast

  !-----------------------------------------------------------------------------
  !> Extended MPI_Ibcast for type real(RNP)

  interface XMPI_Ibcast
    module procedure IbcastX0
    module procedure IbcastX1
    module procedure IbcastX2
  end interface XMPI_Ibcast

  !-----------------------------------------------------------------------------
  !> Extended MPI_Isend for type real(RNP)

  interface XMPI_Isend
    module procedure IsendX0
    module procedure IsendX1
    module procedure IsendX2
    module procedure IsendX3
    module procedure IsendX4
    module procedure IsendX5
  end interface XMPI_Isend

  !-----------------------------------------------------------------------------
  !> Extended MPI_Irecv for type real(RNP)

  interface XMPI_Irecv
    module procedure IrecvX0
    module procedure IrecvX1
    module procedure IrecvX2
    module procedure IrecvX3
    module procedure IrecvX4
    module procedure IrecvX5
  end interface XMPI_Irecv

  !-----------------------------------------------------------------------------
  !> Extended MPI_Reduce for type real(RNP)

  interface XMPI_Reduce
    module procedure ReduceX00
    module procedure ReduceX11
  end interface XMPI_Reduce

  !-----------------------------------------------------------------------------
  !> Extended MPI_Allreduce for type real(RNP)

  interface XMPI_Allreduce
    module procedure AllreduceX00
    module procedure AllreduceX11
  end interface XMPI_Allreduce

contains

!===============================================================================
! MPI_Bcast

!-------------------------------------------------------------------------------
!> Bcast for scalar buffers

subroutine BcastX0(buffer, root, comm)
  real(RNP),      intent(inout) :: buffer !< buffer
  integer,        intent(in)    :: root   !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm   !< communicator

  call MPI_Bcast(buffer, 1, MPI_REAL_RNP, root, comm)

end subroutine BcastX0

!-------------------------------------------------------------------------------
!> Bcast for 1D buffers

subroutine BcastX1(buffer, root, comm)
  real(RNP),      intent(inout) :: buffer(:) !< buffer
  integer,        intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm      !< communicator

  call MPI_Bcast(buffer, size(buffer), MPI_REAL_RNP, root, comm)

end subroutine BcastX1

!-------------------------------------------------------------------------------
!> Bcast for 2D buffers

subroutine BcastX2(buffer, root, comm)
  real(RNP),      intent(inout) :: buffer(:,:) !< buffer
  integer,        intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm        !< communicator

  call MPI_Bcast(buffer, size(buffer), MPI_REAL_RNP, root, comm)

end subroutine BcastX2

!===============================================================================
! MPI_Ibcast

!-------------------------------------------------------------------------------
!> Ibcast for scalar buffers

subroutine IbcastX0(buffer, root, comm, request)
  real(RNP),         intent(inout) :: buffer   !< buffer
  integer,           intent(in)    :: root     !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm     !< communicator
  type(MPI_Request), intent(out)   :: request  !< request

  call MPI_Ibcast(buffer, 1, MPI_REAL_RNP, root, comm, request)

end subroutine IbcastX0

!-------------------------------------------------------------------------------
!> Ibcast for 1D buffers

subroutine IbcastX1(buffer, root, comm, request)
  real(RNP),         intent(inout) :: buffer(:) !< buffer
  integer,           intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm      !< communicator
  type(MPI_Request), intent(out)   :: request   !< request

  call MPI_Ibcast(buffer, size(buffer), MPI_REAL_RNP, root, comm, request)

end subroutine IbcastX1

!-------------------------------------------------------------------------------
!> Ibcast for 2D buffers

subroutine IbcastX2(buffer, root, comm, request)
  real(RNP),         intent(inout) :: buffer(:,:) !< buffer
  integer,           intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm),    intent(in)    :: comm        !< communicator
  type(MPI_Request), intent(out)   :: request     !< request

  call MPI_Ibcast(buffer, size(buffer), MPI_REAL_RNP, root, comm, request)

end subroutine IbcastX2

!===============================================================================
! MPI_Isend

!-------------------------------------------------------------------------------
!> Isend for scalar buffers

subroutine IsendX0(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer   !< send buffer
  integer,           intent(in)  :: dest     !< rank of the receiver
  integer,           intent(in)  :: tag      !< message tag
  type(MPI_Comm),    intent(in)  :: comm     !< communicator
  type(MPI_Request), intent(out) :: request  !< request

  call MPI_Isend(buffer, 1, MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX0

!-------------------------------------------------------------------------------
!> Isend for 1D buffers

subroutine IsendX1(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer(:) !< send buffer
  integer,           intent(in)  :: dest      !< rank of the receiver
  integer,           intent(in)  :: tag       !< message tag
  type(MPI_Comm),    intent(in)  :: comm      !< communicator
  type(MPI_Request), intent(out) :: request   !< request

  call MPI_Isend(buffer, size(buffer), MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX1

!-------------------------------------------------------------------------------
!> Isend for 2D buffers

subroutine IsendX2(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer(:,:) !< send buffer
  integer,           intent(in)  :: dest        !< rank of the receiver
  integer,           intent(in)  :: tag         !< message tag
  type(MPI_Comm),    intent(in)  :: comm        !< communicator
  type(MPI_Request), intent(out) :: request     !< request

  call MPI_Isend(buffer, size(buffer), MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX2

!-------------------------------------------------------------------------------
!> Isend for 3D buffers

subroutine IsendX3(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer(:,:,:) !< send buffer
  integer,           intent(in)  :: dest          !< rank of the receiver
  integer,           intent(in)  :: tag           !< message tag
  type(MPI_Comm),    intent(in)  :: comm          !< communicator
  type(MPI_Request), intent(out) :: request       !< request

  call MPI_Isend(buffer, size(buffer), MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX3

!-------------------------------------------------------------------------------
!> Isend for 4D buffers

subroutine IsendX4(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer(:,:,:,:) !< send buffer
  integer,           intent(in)  :: dest            !< rank of the receiver
  integer,           intent(in)  :: tag             !< message tag
  type(MPI_Comm),    intent(in)  :: comm            !< communicator
  type(MPI_Request), intent(out) :: request         !< request

  call MPI_Isend(buffer, size(buffer), MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX4

!-------------------------------------------------------------------------------
!> Isend for 5D buffers

subroutine IsendX5(buffer, dest, tag, comm, request)
  real(RNP),         intent(in)  :: buffer(:,:,:,:,:) !< send buffer
  integer,           intent(in)  :: dest              !< rank of receiver
  integer,           intent(in)  :: tag               !< message tag
  type(MPI_Comm),    intent(in)  :: comm              !< communicator
  type(MPI_Request), intent(out) :: request           !< request

  call MPI_Isend(buffer, size(buffer), MPI_REAL_RNP, dest, tag, comm, request)

end subroutine IsendX5

!===============================================================================
! MPI_Irecv

!-------------------------------------------------------------------------------
!> Irecv for scalar buffers

subroutine IrecvX0(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer   !< reveive buffer
  integer,           intent(in)  :: source   !< rank of the sender
  integer,           intent(in)  :: tag      !< message tag
  type(MPI_Comm),    intent(in)  :: comm     !< communicator
  type(MPI_Request), intent(out) :: request  !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, 1, MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX0

!-------------------------------------------------------------------------------
!> Irecv for 1D buffers

subroutine IrecvX1(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer(:) !< reveive buffer
  integer,           intent(in)  :: source    !< rank of the sender
  integer,           intent(in)  :: tag       !< message tag
  type(MPI_Comm),    intent(in)  :: comm      !< communicator
  type(MPI_Request), intent(out) :: request   !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, size(buffer), MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX1

!-------------------------------------------------------------------------------
!> Irecv for 2D buffers

subroutine IrecvX2(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer(:,:) !< reveive buffer
  integer,           intent(in)  :: source      !< rank of the sender
  integer,           intent(in)  :: tag         !< message tag
  type(MPI_Comm),    intent(in)  :: comm        !< communicator
  type(MPI_Request), intent(out) :: request     !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, size(buffer), MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX2

!-------------------------------------------------------------------------------
!> Irecv for 3D buffers

subroutine IrecvX3(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer(:,:,:) !< reveive buffer
  integer,           intent(in)  :: source        !< rank of the sender
  integer,           intent(in)  :: tag           !< message tag
  type(MPI_Comm),    intent(in)  :: comm          !< communicator
  type(MPI_Request), intent(out) :: request       !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, size(buffer), MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX3

!-------------------------------------------------------------------------------
!> Irecv for 4D buffers

subroutine IrecvX4(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer(:,:,:,:) !< reveive buffer
  integer,           intent(in)  :: source          !< rank of the sender
  integer,           intent(in)  :: tag             !< message tag
  type(MPI_Comm),    intent(in)  :: comm            !< communicator
  type(MPI_Request), intent(out) :: request         !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, size(buffer), MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX4

!-------------------------------------------------------------------------------
!> Irecv for 5D buffers

subroutine IrecvX5(buffer, source, tag, comm, request)
  real(RNP),         intent(out) :: buffer(:,:,:,:,:) !< reveive buffer
  integer,           intent(in)  :: source            !< rank of the sender
  integer,           intent(in)  :: tag               !< message tag
  type(MPI_Comm),    intent(in)  :: comm              !< communicator
  type(MPI_Request), intent(out) :: request           !< request
  !!! NOT YET SUPPORTED BY PGI !!! asynchronous :: buffer

  call MPI_Irecv(buffer, size(buffer), MPI_REAL_RNP, source, tag, comm, request)

end subroutine IrecvX5

!===============================================================================
! MPI_Reduce

!-------------------------------------------------------------------------------
!> Reduce for scalar/scalar send/receive buffers

subroutine ReduceX00(sendbuf, recvbuf, op, root, comm)
  real(RNP),      intent(in)    :: sendbuf !< send buffer
  real(RNP),      intent(inout) :: recvbuf !< receive buffer
  type(MPI_Op),   intent(in)    :: op      !< reduce operation
  integer,        intent(in)    :: root    !< rank of root process
  type(MPI_Comm), intent(in)    :: comm    !< communicator

  call MPI_Reduce(sendbuf, recvbuf, 1, MPI_REAL_RNP, op, root, comm)

end subroutine ReduceX00

!-------------------------------------------------------------------------------
!> Reduce for 1D/1D send/receive buffers

subroutine ReduceX11(sendbuf, recvbuf, op, root, comm)
  real(RNP),      intent(in)    :: sendbuf(:)             !< send buffer
  real(RNP),      intent(inout) :: recvbuf(size(sendbuf)) !< receive buffer
  type(MPI_Op),   intent(in)    :: op                     !< reduce operation
  integer,        intent(in)    :: root                   !< rank of root
  type(MPI_Comm), intent(in)    :: comm                   !< communicator

  call MPI_Reduce(sendbuf, recvbuf, size(sendbuf), MPI_REAL_RNP, op, root, comm)

end subroutine ReduceX11

!===============================================================================
! MPI_Allreduce

!-------------------------------------------------------------------------------
!> Allreduce for scalar/scalar send/receive buffers

subroutine AllreduceX00(sendbuf, recvbuf, op, comm)
  real(RNP),      intent(in)    :: sendbuf !< send buffer
  real(RNP),      intent(inout) :: recvbuf !< receive buffer
  type(MPI_Op),   intent(in)    :: op      !< reduce operation
  type(MPI_Comm), intent(in)    :: comm    !< communicator

  call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_REAL_RNP, op, comm)

end subroutine AllreduceX00

!-------------------------------------------------------------------------------
!> Allreduce for 1D/1D send/receive buffers

subroutine AllreduceX11(sendbuf, recvbuf, op, comm)
  real(RNP),      intent(in)    :: sendbuf(:)             !< send buffer
  real(RNP),      intent(inout) :: recvbuf(size(sendbuf)) !< receive buffer
  type(MPI_Op),   intent(in)    :: op                     !< reduce operation
  type(MPI_Comm), intent(in)    :: comm                   !< communicator

  call MPI_Allreduce(sendbuf, recvbuf, size(sendbuf), MPI_REAL_RNP, op, comm)

end subroutine AllreduceX11

!===============================================================================

end module XMPI__Real_RNP
