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

!> summary:  Extended MPI Fortran binding for type integer
!> author:   Joerg Stiller
!> date:     2014/11/09, revised 2017/04/10
!===============================================================================

module XMPI__Integer
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
  !> Extended MPI_Bcast for type integer

  interface XMPI_Bcast
    module procedure BcastX0
    module procedure BcastX1
    module procedure BcastX2
    module procedure BcastX3
  end interface XMPI_Bcast

  !-----------------------------------------------------------------------------
  !> Extended MPI_Ibcast for type integer

  interface XMPI_Ibcast
    module procedure IbcastX0
    module procedure IbcastX1
    module procedure IbcastX2
  end interface XMPI_Ibcast

  !-----------------------------------------------------------------------------
  !> Enhanced MPI_Isend for type integer

  interface XMPI_Isend
    module procedure IsendX0
    module procedure IsendX1
    module procedure IsendX2
    module procedure IsendX3
    module procedure IsendX4
  end interface XMPI_Isend

  !-----------------------------------------------------------------------------
  !> Enhanced MPI_Irecv for type integer

  interface XMPI_Irecv
    module procedure IrecvX0
    module procedure IrecvX1
    module procedure IrecvX2
    module procedure IrecvX3
    module procedure IrecvX4
  end interface XMPI_Irecv

  !-----------------------------------------------------------------------------
  !> Extended MPI_Reduce for type integer

  interface XMPI_Reduce
    module procedure ReduceX00
    module procedure ReduceX11
    module procedure ReduceX22
  end interface XMPI_Reduce

  !-----------------------------------------------------------------------------
  !> Extended MPI_Allreduce for type integer

  interface XMPI_Allreduce
    module procedure AllreduceX00
    module procedure AllreduceX11
    module procedure AllreduceX22
  end interface XMPI_Allreduce

contains

!===============================================================================
! MPI_Bcast

!-------------------------------------------------------------------------------
!> Bcast for scalar buffers

subroutine BcastX0(buffer, root, comm)
  integer,        intent(inout) :: buffer !< buffer
  integer,        intent(in)    :: root   !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm   !< communicator

  call MPI_Bcast(buffer, 1, MPI_INTEGER, root, comm)

end subroutine BcastX0

!-------------------------------------------------------------------------------
!> Bcast for 1D buffers

subroutine BcastX1(buffer, root, comm)
  integer,        intent(inout) :: buffer(:) !< buffer
  integer,        intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm      !< communicator

  call MPI_Bcast(buffer, size(buffer), MPI_INTEGER, root, comm)

end subroutine BcastX1

!-------------------------------------------------------------------------------
!> Bcast for 2D buffers

subroutine BcastX2(buffer, root, comm)
  integer,        intent(inout) :: buffer(:,:) !< buffer
  integer,        intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm        !< communicator

  call MPI_Bcast(buffer, size(buffer), MPI_INTEGER, root, comm)

end subroutine BcastX2

!-------------------------------------------------------------------------------
!> Bcast for 3D buffers

subroutine BcastX3(buffer, root, comm)
  integer,        intent(inout) :: buffer(:,:,:) !< buffer
  integer,        intent(in)    :: root          !< rank of broadcast root
  type(MPI_Comm), intent(in)    :: comm          !< communicator

  call MPI_Bcast(buffer, size(buffer), MPI_INTEGER, root, comm)

end subroutine BcastX3

!===============================================================================
! MPI_Ibcast

!-------------------------------------------------------------------------------
!> Ibcast for scalar buffers

subroutine IbcastX0(buffer, root, comm, request)
  integer, asynchronous, intent(inout) :: buffer   !< buffer
  integer,               intent(in)    :: root     !< rank of broadcast root
  type(MPI_Comm),        intent(in)    :: comm     !< communicator
  type(MPI_Request),     intent(out)   :: request  !< request

  call MPI_Ibcast(buffer, 1, MPI_INTEGER, root, comm, request)

end subroutine IbcastX0

!-------------------------------------------------------------------------------
!> Ibcast for 1D buffers

subroutine IbcastX1(buffer, root, comm, request)
  integer, asynchronous, intent(inout) :: buffer(:) !< buffer
  integer,               intent(in)    :: root      !< rank of broadcast root
  type(MPI_Comm),        intent(in)    :: comm      !< communicator
  type(MPI_Request),     intent(out)   :: request   !< request

  call MPI_Ibcast(buffer(1), size(buffer), MPI_INTEGER, root, comm, request)

end subroutine IbcastX1

!-------------------------------------------------------------------------------
!> Ibcast for 2D buffers

subroutine IbcastX2(buffer, root, comm, request)
  integer, asynchronous, intent(inout) :: buffer(:,:) !< buffer
  integer,               intent(in)    :: root        !< rank of broadcast root
  type(MPI_Comm),        intent(in)    :: comm        !< communicator
  type(MPI_Request),     intent(out)   :: request     !< request

  call MPI_Ibcast(buffer(1,1), size(buffer), MPI_INTEGER, root, comm, request)

end subroutine IbcastX2

!===============================================================================
! MPI_Isend

!-------------------------------------------------------------------------------
!> Isend for scalar buffers

subroutine IsendX0(buffer, dest, tag, comm, request)
  integer, asynchronous, intent(in)  :: buffer   !< send buffer
  integer,               intent(in)  :: dest     !< rank of receiver
  integer,               intent(in)  :: tag      !< message tag
  type(MPI_Comm),        intent(in)  :: comm     !< communicator
  type(MPI_Request),     intent(out) :: request  !< request

  call MPI_Isend(buffer, 1, MPI_INTEGER, dest, tag, comm, request)

end subroutine IsendX0

!-------------------------------------------------------------------------------
!> Isend for 1D buffers

subroutine IsendX1(buffer, dest, tag, comm, request)
  integer, asynchronous, intent(in)  :: buffer(:) !< send buffer
  integer,               intent(in)  :: dest      !< rank of receiver
  integer,               intent(in)  :: tag       !< message tag
  type(MPI_Comm),        intent(in)  :: comm      !< communicator
  type(MPI_Request),     intent(out) :: request   !< request

  call MPI_Isend(buffer, size(buffer), MPI_INTEGER, dest, tag, comm, request)

end subroutine IsendX1

!-------------------------------------------------------------------------------
!> Isend for 2D buffers

subroutine IsendX2(buffer, dest, tag, comm, request)
  integer, asynchronous, intent(in)  :: buffer(:,:) !< send buffer
  integer,               intent(in)  :: dest        !< rank of receiver
  integer,               intent(in)  :: tag         !< message tag
  type(MPI_Comm),        intent(in)  :: comm        !< communicator
  type(MPI_Request),     intent(out) :: request     !< request

  call MPI_Isend(buffer, size(buffer), MPI_INTEGER, dest, tag, comm, request)

end subroutine IsendX2

!-------------------------------------------------------------------------------
!> Isend for 3D buffers

subroutine IsendX3(buffer, dest, tag, comm, request)
  integer, asynchronous, intent(in)  :: buffer(:,:,:) !< send buffer
  integer,               intent(in)  :: dest          !< rank of receiver
  integer,               intent(in)  :: tag           !< message tag
  type(MPI_Comm),        intent(in)  :: comm          !< communicator
  type(MPI_Request),     intent(out) :: request       !< request

  call MPI_Isend(buffer, size(buffer), MPI_INTEGER, dest, tag, comm, request)

end subroutine IsendX3

!-------------------------------------------------------------------------------
!> Isend for 4D buffers

subroutine IsendX4(buffer, dest, tag, comm, request)
  integer, asynchronous, intent(in)  :: buffer(:,:,:,:) !< send buffer
  integer,               intent(in)  :: dest            !< rank of receiver
  integer,               intent(in)  :: tag             !< message tag
  type(MPI_Comm),        intent(in)  :: comm            !< communicator
  type(MPI_Request),     intent(out) :: request         !< request

  call MPI_Isend(buffer, size(buffer), MPI_INTEGER, dest, tag, comm, request)

end subroutine IsendX4

!===============================================================================
! MPI_Irecv

!-------------------------------------------------------------------------------
!> Irecv for scalar buffers

subroutine IrecvX0(buffer, source, tag, comm, request)
  integer, asynchronous, intent(out) :: buffer   !< reveive buffer
  integer,               intent(in)  :: source   !< rank of sender
  integer,               intent(in)  :: tag      !< message tag
  type(MPI_Comm),        intent(in)  :: comm     !< communicator
  type(MPI_Request),     intent(out) :: request  !< request

  call MPI_Irecv(buffer, 1, MPI_INTEGER, source, tag, comm, request)

end subroutine IrecvX0

!-------------------------------------------------------------------------------
!> Irecv for 1D buffers

subroutine IrecvX1(buffer, source, tag, comm, request)
  integer, asynchronous, intent(out) :: buffer(:) !< reveive buffer
  integer,               intent(in)  :: source    !< rank of sender
  integer,               intent(in)  :: tag       !< message tag
  type(MPI_Comm),        intent(in)  :: comm      !< communicator
  type(MPI_Request),     intent(out) :: request   !< request

  call MPI_Irecv(buffer, size(buffer), MPI_INTEGER, source, tag, comm, request)

end subroutine IrecvX1

!-------------------------------------------------------------------------------
!> Irecv for 2D buffers

subroutine IrecvX2(buffer, source, tag, comm, request)
  integer, asynchronous, intent(out) :: buffer(:,:) !< reveive buffer
  integer,               intent(in)  :: source      !< rank of sender
  integer,               intent(in)  :: tag         !< message tag
  type(MPI_Comm),        intent(in)  :: comm        !< communicator
  type(MPI_Request),     intent(out) :: request     !< request

  call MPI_Irecv(buffer, size(buffer), MPI_INTEGER, source, tag, comm, request)

end subroutine IrecvX2

!-------------------------------------------------------------------------------
!> Irecv for 3D buffers

subroutine IrecvX3(buffer, source, tag, comm, request)
  integer, asynchronous, intent(out) :: buffer(:,:,:) !< reveive buffer
  integer,               intent(in)  :: source        !< rank of sender
  integer,               intent(in)  :: tag           !< message tag
  type(MPI_Comm),        intent(in)  :: comm          !< communicator
  type(MPI_Request),     intent(out) :: request       !< request

  call MPI_Irecv(buffer, size(buffer), MPI_INTEGER, source, tag, comm, request)

end subroutine IrecvX3

!-------------------------------------------------------------------------------
!> Irecv for 4D buffers

subroutine IrecvX4(buffer, source, tag, comm, request)
  integer, asynchronous, intent(out) :: buffer(:,:,:,:) !< reveive buffer
  integer,               intent(in)  :: source          !< rank of sender
  integer,               intent(in)  :: tag             !< message tag
  type(MPI_Comm),        intent(in)  :: comm            !< communicator
  type(MPI_Request),     intent(out) :: request         !< request

  call MPI_Irecv(buffer, size(buffer), MPI_INTEGER, source, tag, comm, request)

end subroutine IrecvX4

!===============================================================================
! MPI_Reduce

!-------------------------------------------------------------------------------
!> Reduce for scalar/scalar send/receive buffers

subroutine ReduceX00(sendbuf, recvbuf, op, root, comm)
  integer,        intent(in)    :: sendbuf !< send buffer
  integer,        intent(inout) :: recvbuf !< receive buffer
  type(MPI_Op),   intent(in)    :: op      !< reduce operation
  integer,        intent(in)    :: root    !< rank of root process
  type(MPI_Comm), intent(in)    :: comm    !< communicator

  call MPI_Reduce(sendbuf, recvbuf, 1, MPI_INTEGER, op, root, comm)

end subroutine ReduceX00

!-------------------------------------------------------------------------------
!> Reduce for 1D/1D send/receive buffers

subroutine ReduceX11(sendbuf, recvbuf, op, root, comm)
  integer,        intent(in)    :: sendbuf(:)             !< send buffer
  integer,        intent(inout) :: recvbuf(size(sendbuf)) !< receive buffer
  type(MPI_Op),   intent(in)    :: op                     !< reduce operation
  integer,        intent(in)    :: root                   !< rank of root
  type(MPI_Comm), intent(in)    :: comm                   !< communicator

  call MPI_Reduce(sendbuf, recvbuf, size(sendbuf), MPI_INTEGER, op, root, comm)

end subroutine ReduceX11

!-------------------------------------------------------------------------------
!> Reduce for 2D/2D send/receive buffers

subroutine ReduceX22(sendbuf, recvbuf, op, root, comm)
  integer,        intent(in)    :: sendbuf(:,:) !< send buffer
  integer,        intent(inout) :: recvbuf(:,:) !< receive buffer
  type(MPI_Op),   intent(in)    :: op           !< reduce operation
  integer,        intent(in)    :: root         !< rank of root
  type(MPI_Comm), intent(in)    :: comm         !< communicator

  call MPI_Reduce(sendbuf, recvbuf, size(sendbuf), MPI_INTEGER, op, root, comm)

end subroutine ReduceX22

!===============================================================================
! MPI_Allreduce

!-------------------------------------------------------------------------------
!> Allreduce for scalar/scalar send/receive buffers

subroutine AllreduceX00(sendbuf, recvbuf, op, comm)
  integer,        intent(in)    :: sendbuf !< send buffer
  integer,        intent(inout) :: recvbuf !< receive buffer
  type(MPI_Op),   intent(in)    :: op      !< reduce operation
  type(MPI_Comm), intent(in)    :: comm    !< communicator

  call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_INTEGER, op, comm)

end subroutine AllreduceX00

!-------------------------------------------------------------------------------
!> Allreduce for 1D/1D send/receive buffers

subroutine AllreduceX11(sendbuf, recvbuf, op, comm)
  integer,        intent(in)    :: sendbuf(:)             !< send buffer
  integer,        intent(inout) :: recvbuf(size(sendbuf)) !< receive buffer
  type(MPI_Op),   intent(in)    :: op                     !< reduce operation
  type(MPI_Comm), intent(in)    :: comm                   !< communicator

  call MPI_Allreduce(sendbuf, recvbuf, size(sendbuf), MPI_INTEGER, op, comm)

end subroutine AllreduceX11

!-------------------------------------------------------------------------------
!> Allreduce for 2D/2D send/receive buffers

subroutine AllreduceX22(sendbuf, recvbuf, op, comm)
  integer,        intent(in)    :: sendbuf(:,:) !< send buffer
  integer,        intent(inout) :: recvbuf(:,:) !< receive buffer
  type(MPI_Op),   intent(in)    :: op           !< reduce operation
  type(MPI_Comm), intent(in)    :: comm         !< communicator

  call MPI_Allreduce(sendbuf, recvbuf, size(sendbuf), MPI_INTEGER, op, comm)

end subroutine AllreduceX22

!===============================================================================

end module XMPI__Integer
