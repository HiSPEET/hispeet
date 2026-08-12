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

!> summary:  Type and methods for transferring linked face data
!> author:   Joerg Stiller
!> date:     2017/04/30
!>
!> @note
!> This module is obsolet and should not be used
!===============================================================================

module CART__Face_Transfer_Buffer

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use XMPI
  use CART__Mesh_Partition

  implicit none
  private

  public :: FaceTransferBuffer

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked face data and metadata

  type FaceTransferData
    real(RNP),         allocatable :: buf(:)     !< message buffer
    type(MPI_Request), allocatable :: request(:) !< requests
  contains
    final :: Delete_FaceTransferData
  end type FaceTransferData

  !-----------------------------------------------------------------------------
  !> Type for transferring linked face data
  !>
  !> The type supports two flavors of face variables
  !>
  !>    *  single variable of the shape `v(np,np,nf)`
  !>    *  variable arrays of the shape `v(np,np,nf,nc)`
  !>
  !> where `v` is of type `real(RNP)` and
  !>
  !>    *  `nf` is the number of local mesh faces
  !>    *  `nc` is the number of components (variables)
  !>
  !> Use with single thread (no OpenMP):
  !>
  !>        type(FaceTransferBuffer), asynchronous :: face_buf
  !>        ...
  !>        ! create and fill buffer, start transfer
  !>        call face_buf % Transfer(mesh, v, tag)
  !>        ...
  !>        ! possibly perform some computations to hide communication costs
  !>        ...
  !>        ! merge received data and wait transfer to finish
  !>        call face_buf % Merge(mesh, v, alpha, beta)
  !>        call face_buf % Finish()
  !>
  !> Use with OpenMP:
  !>
  !>    *  the buffer must be declared `save` to be shared among the threads
  !>    *  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>
  !>         type(FaceTransferBuffer), asynchronous, allocatable, save :: face_buf
  !>
  !>         !$omp single
  !>         allocate(face_buf)
  !>         !$omp end single
  !>         ...
  !>         call face_buf % Transfer(mesh, v, tag)
  !>         ...
  !>         call face_buf % Merge(mesh, v, alpha, beta)
  !>         call face_buf % Finish()
  !>         ...
  !>         !$omp barrier
  !>         !$omp single
  !>         deallocate(face_buf)
  !>         !$omp end single

  type FaceTransferBuffer
!    private

    integer :: np = 0                  !< number of points per face
    integer :: nc = 0                  !< number of components (variables)
    integer, allocatable   :: start(:) !< message start addresses
    integer, allocatable   :: len(:)   !< message lengths
    type(FaceTransferData) :: send     !< send buffer and metadata
    type(FaceTransferData) :: recv     !< receive buffer and metadata

  contains

    generic,   public  :: New => New_TransferBuffer_S, New_TransferBuffer_A
    procedure, private :: New_TransferBuffer_S
    procedure, private :: New_TransferBuffer_A

    generic,   public  :: Transfer => Transfer_S, Transfer_A
    procedure, private :: Transfer_S
    procedure, private :: Transfer_A

    generic :: Merge => Merge_S, Merge_A
    procedure, private :: Merge_S
    procedure, private :: Merge_A

    procedure :: Finish   => FinishTransfer

    final :: Delete_TransferBuffer

  end type FaceTransferBuffer

contains

!===============================================================================
! FaceTransferData: type bound procedures

!-------------------------------------------------------------------------------
!> Delete face transfer data

subroutine Delete_FaceTransferData(this)
  type(FaceTransferData), intent(inout) :: this  !< buffer

  if (allocated(this % buf    )) deallocate(this % buf    )
  if (allocated(this % request)) deallocate(this % request)

end subroutine Delete_FaceTransferData

!===============================================================================
! FaceTransferBuffer: type bound procedures

!-------------------------------------------------------------------------------
!> Create a new face transfer buffer from mesh and given variable dimensions

subroutine New_TransferBuffer_X(this, mesh, np, nc)
  class(FaceTransferBuffer), intent(inout) :: this     !< buffer
  type(MeshPartition),       intent(in)    :: mesh     !< mesh partition
  integer,                   intent(in)    :: np       !< num points/face
  integer,                   intent(in)    :: nc       !< num components

  integer :: lb, nl
  integer :: i

  !$omp single

  if (allocated(this % start)) deallocate(this % start)
  if (allocated(this % len  )) deallocate(this % len  )

  this%np = np
  this%nc = nc

  nl = size(mesh%link)

  ! start adresses and lengths .................................................

  allocate( this%start(nl) )
  allocate( this%len(nl) )

  associate(m => this%start, l => this%len)
    lb = 0
    do i = 1, nl
      if (i == 1) then
        m(i) = 1
      else
        m(i) = m(i-1) + l(i-1)
      end if
      l(i) = np * nc * mesh%link(i)%nf
      lb = lb + l(i)
    end do
  end associate

  ! buffers ....................................................................

  allocate(this % send % buf( lb ))
  allocate(this % recv % buf( lb ))

  ! requests ..................................................................

  allocate( this % send % request(nl), source = MPI_REQUEST_NULL )
  allocate( this % recv % request(nl), source = MPI_REQUEST_NULL )

  !$omp end single

end subroutine New_TransferBuffer_X

!-------------------------------------------------------------------------------
!> Create a new face transfer buffer for a single variable

subroutine New_TransferBuffer_S(this, mesh, v)
  class(FaceTransferBuffer), intent(inout) :: this      !< buffer
  type(MeshPartition),       intent(in)    :: mesh      !< mesh partition
  real(RNP),                 intent(in)    :: v(:,:,:)  !< face variable

  integer :: np

  np = size(v,1) * size(v,2)

  call New_TransferBuffer_X(this, mesh, np, nc=1)

end subroutine New_TransferBuffer_S

!-------------------------------------------------------------------------------
!> Create a new face transfer buffer for an array of variables

subroutine New_TransferBuffer_A(this, mesh, v)
  class(FaceTransferBuffer), intent(inout) :: this        !< buffer
  type(MeshPartition),       intent(in)    :: mesh        !< mesh partition
  real(RNP),                 intent(in)    :: v(:,:,:,:)  !< face variable

  integer :: np, nc

  np = size(v,1) * size(v,2)
  nc = size(v,4)

  call New_TransferBuffer_X(this, mesh, np, nc)

end subroutine New_TransferBuffer_A

!-------------------------------------------------------------------------------
!> Extract and transfer buffer -- eXplicit shape version

subroutine Transfer_X(this, mesh, v, tag)
  class(FaceTransferBuffer), asynchronous, intent(inout) :: this  !< buffer
  type(MeshPartition), intent(in) :: mesh  !< mesh partition
  real(RNP),           intent(in) :: v     !< face variable
  integer,             intent(in) :: tag   !< message tag
  dimension :: v(this%np, mesh%nf,this%nc)

  integer :: i, l, m, b1, part

  associate(send => this%send, recv => this%recv, comm => mesh%comm)

    if (size(send%buf) < 1) return

    ! extract send buffer ......................................................

    b1 = 1  ! first entry in send%buf

    do l = 1, size(mesh%link)

      if (mesh%link(l)%part /= mesh%part) then

        call CopyToBuffer(  np   = this%np            &
                         ,  nc   = this%nc            &
                         ,  nf   = mesh%link(l)%nf    &
                         ,  nm   = mesh%nf            &
                         ,  face = mesh%link(l)%face  &
                         ,  v    = v                  &
                         ,  vb   = send%buf(b1:)      )

      else ! local link: copy coupled face data into send buffer

        call CopyToBuffer(  np   = this%np                    &
                         ,  nc   = this%nc                    &
                         ,  nf   = mesh%link(l)%nf            &
                         ,  nm   = mesh%nf                    &
                         ,  face = mesh%link(l)%coupled_face  &
                         ,  v    = v                          &
                         ,  vb   = send%buf(b1:)              )

      end if

      b1 = b1 + this%np * mesh%link(l)%nf * this%nc

    end do

    ! start send/receive .......................................................

    !$omp barrier
    !$omp master

    do i = 1, size(mesh%link)

      part = mesh % link(i) % part

      m = this % start(i)
      l = this % len(i)
      if (l < 1) cycle

      if (part /= mesh%part) then

        call MPI_Isend( send%buf(m:), l, MPI_REAL_RNP, part, tag, comm, &
                        send%request(i)                                  )

        call MPI_Irecv( recv%buf(m:), l, MPI_REAL_RNP, part, tag, comm, &
                        recv%request(i)                                 )

      else
        recv%buf(m:m+l-1) = send%buf(m:m+l-1)
      end if

    end do

    !$omp end master

  end associate

contains

  subroutine CopyToBuffer(np, nc, nf, nm, face, v, vb)
    integer,   intent(in)    :: np            !< number of points/face
    integer,   intent(in)    :: nc            !< number of components
    integer,   intent(in)    :: nf            !< number of linked faces
    integer,   intent(in)    :: nm            !< number of mesh faces
    integer,   intent(in)    :: face(nf)      !< list of linked faces
    real(RNP), intent(in)    :: v(np,nm,nc)   !< face variable
    real(RNP), intent(inout) :: vb(np*nf*nc)  !< buffer

    integer :: i, j, k

    !$omp do collapse(2) private(i,j,k)
    !$acc parallel loop collapse(3) present(v) copyin(face) copyout(vb)
    do k = 1, nc
    do j = 1, nf
    do i = 1, np
      vb(i + np*(j-1 + nf*(k-1))) = v(i, face(j), k)
    end do
    end do
    end do

  end subroutine CopyToBuffer

end subroutine Transfer_X

!-------------------------------------------------------------------------------
!> Extract and transfer real buffer -- single variable
!>
!> The variable must be dimensioned `v(n1,n2,mesh%nf)` with `n1*n2 = this%np`

subroutine Transfer_S(this, mesh, v, tag)
  class(FaceTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition), intent(in) :: mesh     !< mesh partition
  real(RNP),           intent(in) :: v(:,:,:) !< face variable
  integer,             intent(in) :: tag      !< message tag

  if (this%np /= size(v,1)*size(v,2) .or. this%nc /= 1) then
    call this % New(mesh, v)
  end if

  call Transfer_X(this, mesh, v, tag)

end subroutine Transfer_S

!-------------------------------------------------------------------------------
!> Extract and transfer real buffer -- array of variables
!>
!> The variable must be dimensioned `v(n1,n2,mesh%nf,this%nc)`
!> with `n1*n2 = this%np`

subroutine Transfer_A(this, mesh, v, tag)
  class(FaceTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition), intent(in) :: mesh       !< mesh partition
  real(RNP),           intent(in) :: v(:,:,:,:) !< face variable
  integer,             intent(in) :: tag        !< message tag

  if (this%np /= size(v,1)*size(v,2) .or. this%nc /= size(v,4)) then
    call this % New(mesh, v)
  end if

  call Transfer_X(this, mesh, v, tag)

end subroutine Transfer_A

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge buffer into face variable
!>
!> Denoting the buffer with vb, the following operation will be executed:
!>
!>      v  =  alpha * v  +  beta * vb

subroutine Merge_X(this, mesh, v, alpha, beta)

  ! arguments ..................................................................

  class(FaceTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition), intent(in)    :: mesh  !< mesh partition
  real(RNP),           intent(inout) :: v     !< face variable
  real(RNP), optional, intent(in)    :: alpha !< coeff of v  [1]
  real(RNP), optional, intent(in)    :: beta  !< coeff of vb [1]

  dimension :: v(this%np, mesh%nf, this%nc)

  ! internal data ..............................................................

  type(MPI_Status), allocatable :: status(:)

  integer :: l, nr, b1
  real(RNP) :: a, b

  ! initialization .............................................................

  if (size(this%recv%buf) < 1) return

  if (present(alpha)) then
    a = alpha
  else
    a = 1
  end if

  if (present(beta)) then
    b = beta
  else
    b = 1
  end if

  ! wait for receive to complete ...............................................

  !$omp master
  nr = size(this%recv%request)
  allocate(status(nr))
  call MPI_Waitall(nr, this%recv%request, status)
  !$omp end master
  !$omp barrier

  ! merge buffer ...............................................................

  b1 = 1  ! first entry in recv%buf

  do l = 1, size(mesh%link)

    call MergeBuffer(  np   = this%np            &
                    ,  nc   = this%nc            &
                    ,  nf   = mesh%link(l)%nf    &
                    ,  nm   = mesh%nf            &
                    ,  face = mesh%link(l)%face  &
                    ,  a    = a                  &
                    ,  b    = b                  &
                    ,  v    = v                  &
                    ,  vb   = this%recv%buf(b1:) )

    b1 = b1 + this%np * mesh%link(l)%nf * this%nc

  end do

contains

  subroutine MergeBuffer(np, nc, nf, nm, face, a, b, v, vb)
    integer,   intent(in)    :: np            !< number of points/face
    integer,   intent(in)    :: nc            !< number of components
    integer,   intent(in)    :: nf            !< number of linked faces
    integer,   intent(in)    :: nm            !< number of mesh faces
    integer,   intent(in)    :: face(nf)      !< list of linked faces
    real(RNP), intent(in)    :: a             !< coefficient of v
    real(RNP), intent(in)    :: b             !< coefficient of vb
    real(RNP), intent(inout) :: v(np,nm,nc)   !< face variable
    real(RNP), intent(in)    :: vb(np*nf*nc)  !< buffer

    integer :: i, j, k

    if (a /= ZERO) then
      !$omp do collapse(2) private(i,j,k)
      do k = 1, nc
      do j = 1, nf
      do i = 1, np
        v(i, face(j), k)  =  a * v(i, face(j), k)  &
                          +  b * vb(i + np*(j-1 + nf*(k-1)))
      end do
      end do
      end do
    else
      !$omp do collapse(2) private(i,j,k)
      do k = 1, nc
      do j = 1, nf
      do i = 1, np
        v(i, face(j), k)  =  b * vb(i + np*(j-1 + nf*(k-1)))
      end do
      end do
      end do
    end if

  end subroutine MergeBuffer

end subroutine Merge_X

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge real buffer -- single variable
!>
!> The variable must be dimensioned `v(n1,n2,mesh%nf)` with `n1*n2 = this%np`

subroutine Merge_S(this, mesh, v, alpha, beta)
  class(FaceTransferBuffer), intent(inout) :: this      !< buffer
  type(MeshPartition),       intent(in)    :: mesh      !< mesh partition
  real(RNP),                 intent(inout) :: v(:,:,:)  !< face variable
  real(RNP),       optional, intent(in)    :: alpha     !< coeff of v  [1]
  real(RNP),       optional, intent(in)    :: beta      !< coeff of vb [1]

  call Merge_X(this, mesh, v, alpha, beta)

end subroutine Merge_S

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge real buffer -- array of variables
!>
!> The variable must be dimensioned `v(n1,n2,mesh%nf,this%nc)`
!> with `n1*n2 = this%np`

subroutine Merge_A(this, mesh, v, alpha, beta)
  class(FaceTransferBuffer), intent(inout) :: this        !< buffer
  type(MeshPartition),       intent(in)    :: mesh        !< mesh partition
  real(RNP),                 intent(inout) :: v(:,:,:,:)  !< face variable
  real(RNP),       optional, intent(in)    :: alpha       !< coeff of v  [1]
  real(RNP),       optional, intent(in)    :: beta        !< coeff of vb [1]

  call Merge_X(this, mesh, v, alpha, beta)

end subroutine Merge_A

!-------------------------------------------------------------------------------
!> Executes MPI_Waitall to complete all transfers associated with the buffer

subroutine FinishTransfer(this)
  class(FaceTransferBuffer), intent(inout) :: this  !< buffer

  type(MPI_Status), allocatable :: status(:)
  integer :: n

  !$omp master

  n = size(this%send%request)

  if (n > 0) then
    allocate(status(n))
    call MPI_Waitall(n, this%send%request, status)
    call MPI_Waitall(n, this%recv%request, status)
  end if

  !$omp end master
  !$omp barrier

end subroutine FinishTransfer

!-------------------------------------------------------------------------------
!> Delete a face transfer buffer

subroutine Delete_TransferBuffer(this)
  type(FaceTransferBuffer), intent(inout) :: this  !< buffer

  this % np = 0
  this % nc = 0

  if (allocated(this % start)) deallocate(this % start)
  if (allocated(this % len  )) deallocate(this % len  )

end subroutine Delete_TransferBuffer

!===============================================================================

end module CART__Face_Transfer_Buffer
