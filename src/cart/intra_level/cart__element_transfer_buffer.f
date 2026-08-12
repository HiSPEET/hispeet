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

!> summary:  Type and methods for transferring linked element data
!> author:   Joerg Stiller
!> date:     2017/09/05
!===============================================================================

module CART__Element_Transfer_Buffer

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO
  use Array_Assignments, only: ScaleArray
  use Execution_Control, only: Error
  use XMPI
  use CART__Mesh_Link
  use CART__Mesh_Partition

  implicit none
  private

  public :: ElementTransferBuffer

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked element data and metadata

  type ElementTransferData
    integer,           allocatable :: start(:)   !< message start addresses
    integer,           allocatable :: len(:)     !< message lengths
    integer,           allocatable :: node(:)    !< mesh point IDs
    real(RNP),         allocatable :: buf(:)     !< message buffer
    type(MPI_Request), allocatable :: request(:) !< MPI requests
  contains
    final :: Delete_ElementTransferData
  end type ElementTransferData

  !-----------------------------------------------------------------------------
  !> Type for transferring linked element data
  !>
  !> The data of linked remote elements is stored in "ghost" elements and,
  !> hence, can be transferred in two directions:
  !>
  !>   *  from master elements to their ghosts in remote partitions, and
  !>   *  from ghosts to the master elements.
  !>
  !> Depending on the transfer direction, the received data is merged into
  !> the ghost or master element data, respectively.
  !> For realizing anisotropic, element-overlapping operations, the number
  !> of element "points", `np`, can vary in each coordinate direction.
  !> Moreover, the implementation supports full and partial transfers.
  !> In the latter case, only `nl` point layers are transferred. The number
  !> of layers can be chosen independently for each coordinate direction.
  !> If `nl` is not specified or identical to zero, a full copy is performed.
  !>
  !> The type supports two flavors of mesh variables
  !>
  !>    *  single variable of the shape `v(np(1),np(2),np(3),ne+ng)`
  !>    *  variable arrays of the shape `v(np(1),np(2),np(3),ne+ng,nc)`
  !>
  !> where `v` is of type `real(RNP)` and
  !>
  !>    *  `np(1:3)` is the number of element or region points per direction
  !>    *  `ne` is the number of local elements, i.e. `mesh%ne`
  !>    *  `ng` is the number of ghost elements, i.e. `mesh%ng`
  !>    *  `nc` is the number of components (variables)
  !>
  !> It is recommended to pass the variable arguments always in the shape
  !> given above, though the ghost entries are not used (and can be omitted)
  !> in a few cases.
  !>
  !> Use with single thread (no OpenMP):
  !>
  !>    *  example for master to ghost transfer
  !>
  !>         type(ElementTransferBuffer), asynchronous :: buf
  !>         ...
  !>         ! create and fill buffer, start transfer
  !>         call buf % New(mesh, np, nc, nl)
  !>         call buf % ToGhost_Transfer(mesh, v, tag)
  !>         ...
  !>         ! possibly perform some computations to hide communication costs
  !>         ...
  !>         ! merge received master data into ghosts and finish transfer
  !>         call buf % ToGhost_Merge(v, alpha, beta)
  !>         call buf % ToGhost_Finish()
  !>
  !> Use with OpenMP:
  !>
  !>    *  the buffer must be declared `save` to be shared among the threads
  !>    *  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>    *  example for ghost to master transfer
  !>
  !>         type(ElementTransferBuffer), asynchronous, allocatable, save :: buf
  !>         ...
  !>         !$omp single
  !>         allocate(buf)
  !>         !$omp end single
  !>         ...
  !>         call buf % New(mesh, np, nc, nl)
  !>         call buf % ToMaster_Transfer(mesh, v, tag)
  !>         ...
  !>         call buf % ToMaster_Merge(v, alpha, beta)
  !>         call buf % ToMaster_Finish()
  !>         ...
  !>         !$omp barrier
  !>         !$omp single
  !>         deallocate(buf)
  !>         !$omp end single

  type ElementTransferBuffer

    integer :: np(3) = 0  !< points per direction
    integer :: ne    = 0  !< number of local (master) elements
    integer :: ng    = 0  !< number of ghost elements
    integer :: nc    = 0  !< number of components (variables)
    integer :: nl(3) = 0  !< number of layers, all(nl < 1) implies full transfer

    type(ElementTransferData) :: master  !< master element data
    type(ElementTransferData) :: ghost   !< ghost element data

  contains

    ! public procedures ........................................................

    generic,   public  :: New => New_TransferBuffer_X, &
                                 New_TransferBuffer_S, &
                                 New_TransferBuffer_A

    generic,   public  :: ToGhost_Transfer  => ToGhost_Transfer_S,  &
                                               ToGhost_Transfer_A

    generic,   public  :: ToGhost_Merge     => ToGhost_Merge_S,     &
                                               ToGhost_Merge_A
    procedure, public  :: ToGhost_Finish

    generic,   public  :: ToMaster_Transfer => ToMaster_Transfer_S, &
                                               ToMaster_Transfer_A

    generic,   public  :: ToMaster_Merge    => ToMaster_Merge_S,    &
                                               ToMaster_Merge_A
    procedure, public  :: ToMaster_Finish

    ! private procedures .......................................................

    procedure, private :: New_TransferBuffer_X
    procedure, private :: New_TransferBuffer_S
    procedure, private :: New_TransferBuffer_A

    procedure, private :: ToGhost_Transfer_S
    procedure, private :: ToGhost_Transfer_A

    procedure, private :: ToGhost_Merge_S
    procedure, private :: ToGhost_Merge_A

    procedure, private :: ToMaster_Transfer_S
    procedure, private :: ToMaster_Transfer_A

    procedure, private :: ToMaster_Merge_S
    procedure, private :: ToMaster_Merge_A

    ! finalization .............................................................

    final :: Delete_TransferBuffer

  end type ElementTransferBuffer

contains

!===============================================================================
! ElementTransferData: type bound procedures

!-------------------------------------------------------------------------------
!> Delete element transfer data

subroutine Delete_ElementTransferData(this)
  type(ElementTransferData), intent(inout) :: this  !< buffer

  if (allocated(this % start  )) deallocate(this % start  )
  if (allocated(this % len    )) deallocate(this % len    )
  if (allocated(this % node   )) deallocate(this % node   )
  if (allocated(this % buf    )) deallocate(this % buf    )
  if (allocated(this % request)) deallocate(this % request)

end subroutine Delete_ElementTransferData

!===============================================================================
! ElementTransferBuffer: type bound procedures

!-------------------------------------------------------------------------------
!> Create a new element transfer buffer from mesh and given dimensions

subroutine New_TransferBuffer_X(this, mesh, np, nc, nl)
  class(ElementTransferBuffer), intent(inout) :: this  !< buffer
  type(MeshPartition),          intent(in)    :: mesh  !< mesh partition
  integer,                      intent(in)    :: np(3) !< points per direction
  integer,                      intent(in)    :: nc    !< number of components
  integer,            optional, intent(in)    :: nl(3) !< number of layers

  logical, allocatable :: mask(:,:,:)
  integer, allocatable :: node(:)

  integer :: lm, lg, np1, np12, np123
  integer :: e, i, j, k, l, m, n
  logical :: partial

  !$omp single

  ! setup ......................................................................

  this % np = np
  this % ne = mesh % ne
  this % ng = mesh % ng
  this % nc = nc
  if (present(nl)) then
    this % nl = nl
  else
    this % nl = 0
  end if
  partial = any(this%nl > 0)

  lm = sum(mesh % link % nm)
  lg = sum(mesh % link % ng)

  np1   = np(1)
  np12  = np(1) * np(2)
  np123 = np(1) * np(2) * np(3)

  allocate( mask(np(1), np(2), np(3)), source = .true. )
  allocate( node(np123 * max(lm, lg)) )

  ! master element data ........................................................

  allocate(this % master % start( size(mesh%link) ))
  allocate(this % master % len  ( size(mesh%link) ))

  associate(start => this%master%start, len => this%master%len)

    n = 0
    do l = 1, size(mesh%link)

      associate(element => mesh%link(l)%master)
        do m = 1, size(element)
          if (any(this%nl > 0)) then
            call GeneratePointMask(element(m), this%np, this%nl, mask)
          end if
          e = element(m) % id
          do k = 1, np(3)
          do j = 1, np(2)
          do i = 1, np(1)
            if (mask(i,j,k)) then
              n = n + 1
              ! lexical mesh-point index
              node(n) = i + np1 * (j-1) + np12 * (k-1) + np123 * (e-1)
            end if
          end do
          end do
          end do
        end do
      end associate

      if (l == 1) then
        start (l) = 1
        len   (l) = n
      else
        start (l) = start(l-1) + len(l-1)
        len   (l) = n + 1      - start(l)
      end if

    end do

    allocate(this % master % node(n), source = node(1:n))
    allocate(this % master % buf (n * nc))

    allocate(this % master % request(size(mesh%link)), source = MPI_REQUEST_NULL)

  end associate

  ! ghost element data .........................................................

  allocate(this % ghost % start( size(mesh%link) ))
  allocate(this % ghost % len  ( size(mesh%link) ))

  associate(start => this%ghost%start, len => this%ghost%len)

    n = 0
    do l = 1, size(mesh%link)

      associate(element => mesh%link(l)%ghost)
        do m = 1, size(element)
          if (partial) then
            call GeneratePointMask(element(m), np, nl, mask)
          end if
          e = element(m) % id
          do k = 1, np(3)
          do j = 1, np(2)
          do i = 1, np(1)
            if (mask(i,j,k)) then
              n = n + 1
              ! lexical mesh-point index
              node(n) = i + np1 * (j-1) + np12 * (k-1) + np123 * (e-1)
            end if
          end do
          end do
          end do
        end do
      end associate

      if (l == 1) then
        start (l) = 1
        len   (l) = n
      else
        start (l) = start(l-1) + len(l-1)
        len   (l) = n + 1      - start(l)
      end if

    end do

    allocate(this % ghost % node(n), source = node(1:n))
    allocate(this % ghost % buf (n * nc))

    allocate(this % ghost % request(size(mesh%link)), source = MPI_REQUEST_NULL)

  end associate

  !$omp end single

end subroutine New_TransferBuffer_X

!-------------------------------------------------------------------------------
!> Generates a mask of points subjected to transfer operations

pure subroutine GeneratePointMask(link, np, nl, mask)
  class(ElementLink), intent(in)  :: link        !< element link
  integer,            intent(in)  :: np(3)       !< points per direction
  integer,            intent(in)  :: nl(3)       !< point layers per direction
  logical,            intent(out) :: mask(:,:,:) !< mask of linked points

  integer :: l1, l2, l3, r1, r2, r3

  l1 = min(nl(1), np(1))
  l2 = min(nl(2), np(2))
  l3 = min(nl(3), np(3))

  r1 = max(np(1) - nl(1) + 1, 1)
  r2 = max(np(2) - nl(2) + 1, 1)
  r3 = max(np(3) - nl(3) + 1, 1)

  mask = .false.

  ! faces
  if (link % face(1) > 0) mask(  :l1,   :  ,   :  ) = .true.
  if (link % face(2) > 0) mask(r1:  ,   :  ,   :  ) = .true.
  if (link % face(3) > 0) mask(  :  ,   :l2,   :  ) = .true.
  if (link % face(4) > 0) mask(  :  , r2:  ,   :  ) = .true.
  if (link % face(5) > 0) mask(  :  ,   :  ,   :l3) = .true.
  if (link % face(6) > 0) mask(  :  ,   :  , r3:  ) = .true.

  ! x1-edges
  if (link % edge( 1) > 0) mask( : ,   :l2,   :l3) = .true.
  if (link % edge( 2) > 0) mask( : , r2:  ,   :l3) = .true.
  if (link % edge( 3) > 0) mask( : ,   :l2, r3:  ) = .true.
  if (link % edge( 4) > 0) mask( : , r2:  , r3:  ) = .true.

  ! x2-edges
  if (link % edge( 5) > 0) mask(  :l1,  : ,   :l3) = .true.
  if (link % edge( 6) > 0) mask(r1:  ,  : ,   :l3) = .true.
  if (link % edge( 7) > 0) mask(  :l1,  : , r3:  ) = .true.
  if (link % edge( 8) > 0) mask(r1:  ,  : , r3:  ) = .true.

  ! x3-edges
  if (link % edge( 9) > 0) mask(  :l1,   :l2,  : ) = .true.
  if (link % edge(10) > 0) mask(r1:  ,   :l2,  : ) = .true.
  if (link % edge(11) > 0) mask(  :l1, r2:  ,  : ) = .true.
  if (link % edge(12) > 0) mask(r1:  , r2:  ,  : ) = .true.

  ! vertices
  if (link % vertex( 1) > 0) mask(  :l1,   :l2,   :l3) = .true.
  if (link % vertex( 2) > 0) mask(r1:  ,   :l2,   :l3) = .true.
  if (link % vertex( 3) > 0) mask(  :l1, r2:  ,   :l3) = .true.
  if (link % vertex( 4) > 0) mask(r1:  , r2:  ,   :l3) = .true.
  if (link % vertex( 5) > 0) mask(  :l1,   :l2, r3:  ) = .true.
  if (link % vertex( 6) > 0) mask(r1:  ,   :l2, r3:  ) = .true.
  if (link % vertex( 7) > 0) mask(  :l1, r2:  , r3:  ) = .true.
  if (link % vertex( 8) > 0) mask(r1:  , r2:  , r3:  ) = .true.

end subroutine GeneratePointMask

!-------------------------------------------------------------------------------
!> Create a new element transfer buffer for a single variable

subroutine New_TransferBuffer_S(this, mesh, v, nl)
  class(ElementTransferBuffer),  intent(inout) :: this   !< buffer
  type(MeshPartition),           intent(in)    :: mesh   !< mesh partition
  real(RNP), dimension(:,:,:,:), intent(in)    :: v      !< mesh variable
  integer,             optional, intent(in)    :: nl(3)  !< number of layers

  integer :: np(3), nc

  np(1) = size(v,1)
  np(2) = size(v,2)
  np(3) = size(v,3)
  nc    = 1

  call New_TransferBuffer_X(this, mesh, np, nc, nl)

end subroutine New_TransferBuffer_S

!-------------------------------------------------------------------------------
!> Create a new element transfer buffer for an array of variables

subroutine New_TransferBuffer_A(this, mesh, v, nl)
  class(ElementTransferBuffer),    intent(inout) :: this   !< buffer
  type(MeshPartition),             intent(in)    :: mesh   !< mesh partition
  real(RNP), dimension(:,:,:,:,:), intent(in)    :: v      !< mesh variable
  integer,               optional, intent(in)    :: nl(3)  !< number of layers

  integer :: np(3), nc

  np(1) = size(v,1)
  np(2) = size(v,2)
  np(3) = size(v,3)
  nc    = size(v,5)

  call New_TransferBuffer_X(this, mesh, np, nc, nl)

end subroutine New_TransferBuffer_A

!===============================================================================

!-------------------------------------------------------------------------------
!> Send master data to ghosts and receive own ghost data -- eXplicit

subroutine ToGhost_Transfer_X(this, mesh, ne_t, v, tag)

  ! arguments ..................................................................

  class(ElementTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition),          intent(in)    :: mesh  !< mesh partition
  integer,                      intent(in)    :: ne_t  !< total num elements
  real(RNP),                    intent(in)    :: v     !< mesh variable
  integer,                      intent(in)    :: tag   !< message tag

  dimension :: v(this%np(1), this%np(2), this%np(3), ne_t, this%nc)

  ! internal data ..............................................................

  integer :: dest, source
  integer :: i, l, m

  ! sanity check ..............................................................

  !$omp master
  if (ne_t < this%ne) then
    call Error( 'ToGhost_Transfer_X'             &
              , 'size(v,4) < this%ne'           &
              , 'CART__Element_Transfer_Buffer' )
  end if
  !$omp end master

  associate(master => this%master, ghost => this%ghost, comm => mesh%comm)

    ! copy master data into buffer ............................................

    call CopyToBuffer( nn   = size(master%node)  &
                     , nm   = size(v) / this%nc  &
                     , nc   = this%nc            &
                     , node = master%node        &
                     , v    = v                  &
                     , vb   = master%buf         )

    ! send master data ........................................................

    !$omp master
    associate(buf => master%buf, request => master%request)
      do i = 1, size(mesh%link)
        dest = mesh % link(i) % part
        m = master % start(i)
        l = master % len(i)
        if (l < 1) cycle
        call MPI_Isend(buf(m:), l, MPI_REAL_RNP, dest, tag, comm, request(i))
      end do
    end associate
    !$omp end master

    ! receive master data ......................................................

    !$omp master
    associate(buf => ghost%buf, request => ghost%request)
      do i = 1, size(mesh%link)
        source = mesh % link(i) % part
        m = ghost % start(i)
        l = ghost % len(i)
        if (l < 1) cycle
        call MPI_Irecv(buf(m:), l, MPI_REAL_RNP, source, tag, comm, request(i))
      end do
    end associate
    !$omp end master

  end associate

contains

  subroutine CopyToBuffer(nn, nm, nc, node, v, vb)
    integer,   intent(in)  :: nn
    integer,   intent(in)  :: nm
    integer,   intent(in)  :: nc
    integer,   intent(in)  :: node(nn)
    real(RNP), intent(in)  :: v(nm,nc)
    real(RNP), intent(out) :: vb(nn,nc)

    integer :: i, j

    !$omp do collapse(2) private(i,j)
    !$acc parallel loop collapse(2) present(v) copyin(node) copyout(vb)
    do j = 1, nc
    do i = 1, nn
      vb(i,j) = v(node(i), j)
    end do
    end do

  end subroutine CopyToBuffer

end subroutine ToGhost_Transfer_X

!-------------------------------------------------------------------------------
!> Send master data to ghosts and receive own ghost data -- single

subroutine ToGhost_Transfer_S(this, mesh, v, tag)
  class(ElementTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition),           intent(in) :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
  integer,                       intent(in) :: tag   !< message tag

  call ToGhost_Transfer_X(this, mesh, size(v,4), v, tag)

end subroutine ToGhost_Transfer_S

!-------------------------------------------------------------------------------
!> Send master data to ghosts and receive own ghost data -- array

subroutine ToGhost_Transfer_A(this, mesh, v, tag)
  class(ElementTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition),             intent(in) :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
  integer,                         intent(in) :: tag   !< message tag

  call ToGhost_Transfer_X(this, mesh, size(v,4), v, tag)

end subroutine ToGhost_Transfer_A

!-------------------------------------------------------------------------------
!> Complete receive merge buffer into ghost data -- eXplicit version
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>       v  =  alpha * v  +  beta * vb

subroutine ToGhost_Merge_X(this, v, alpha, beta)

  ! arguments ..................................................................

  class(ElementTransferBuffer), intent(inout) :: this !< buffer
  real(RNP),           intent(inout) :: v      !< mesh variable
  real(RNP), optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP), optional, intent(in)    :: beta   !< coeff of vb [1]

  dimension :: v(this%np(1), this%np(2), this%np(3), this%ne+this%ng, this%nc)

  ! internal data ..............................................................

  type(MPI_Status), allocatable :: status(:)

  real(RNP) :: a, b
  integer   :: n

  ! initialization .............................................................

  if (size(this % ghost % buf) == 0) return

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
  n = size(this % ghost % request)
  allocate(status(n))
  call MPI_Waitall(n, this%ghost%request, status)
  !$omp end master
  !$omp barrier

  ! merge buffer ...............................................................

  call MergeBuffer( nn   = size(this%ghost%node)  &
                  , nm   = size(v) / this%nc      &
                  , nc   = this%nc                &
                  , node = this % ghost % node    &
                  , a    = a                      &
                  , b    = b                      &
                  , vb   = this%ghost%buf         &
                  , v    = v                      )

contains

  subroutine MergeBuffer(nn, nm, nc, node, a, b, vb, v)
    integer,   intent(in)    :: nn
    integer,   intent(in)    :: nm
    integer,   intent(in)    :: nc
    integer,   intent(in)    :: node(nn)
    real(RNP), intent(in)    :: a, b
    real(RNP), intent(in)    :: vb(nn,nc)
    real(RNP), intent(inout) :: v(nm,nc)

    integer :: i, j

    if (a /= ZERO) then
      !$omp do collapse(2) private(i,j)
      !$acc parallel collapse(2) loop present(v) copyin(node,vb)
      do j = 1, nc
      do i = 1, nn
        v(node(i), j) = a * v(node(i), j)  +  b * vb(i,j)
      end do
      end do
    else
      !$omp do collapse(2) private(i,j)
      !$acc parallel collapse(2) loop present(v) copyin(node,vb)
      do j = 1, nc
      do i = 1, nn
        v(node(i), j) = b * vb(i,j)
      end do
      end do
    end if

  end subroutine MergeBuffer

end subroutine ToGhost_Merge_X

!-------------------------------------------------------------------------------
!> Complete receive merge buffer into ghost data -- single variable
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>     v  =  alpha * v  +  beta * vb
!>
!> As this operation affects the ghost elements, the mesh variable must
!> be dimensioned as `v(np(1), np(2), np(3), ne+ng)`, where
!>
!>    *  `np` is the number of points per direction, as in this%np
!>    *  `ne` is the number of master elements
!>    *  `ng` is the number of ghost elements

subroutine ToGhost_Merge_S(this, v, alpha, beta)
  class(ElementTransferBuffer), intent(inout) :: this !< buffer
  real(RNP), dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
  real(RNP),           optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP),           optional, intent(in)    :: beta   !< coeff of vb [1]

  !$omp master
  if (size(v,4) /= this%ne + this%ng) then
    call Error( 'ToGhost_Merge_S'                &
              , 'size(v,4) /= this%ne + this%ng' &
              , 'CART__Element_Transfer_Buffer'  )
  end if
  !$omp end master

  call ToGhost_Merge_X(this, v, alpha, beta)

end subroutine ToGhost_Merge_S

!-------------------------------------------------------------------------------
!> Complete receive merge buffer into ghost data -- array variable
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>      v  =  alpha * v  +  beta * vb
!>
!> As this operation affects the ghost elements, the mesh variable must
!> be dimensioned as `v(np(1), np(2), np(3), ne+ng, nc)`, where
!>
!>    *  `np` is the number of points per direction, as in `this%np`
!>    *  `ne` is the number of master elements
!>    *  `ng` is the number of ghost elements
!>    *  `nc` is the number of components, as in `this%nc`

subroutine ToGhost_Merge_A(this, v, alpha, beta)
  class(ElementTransferBuffer), intent(inout) :: this !< buffer
  real(RNP), dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
  real(RNP),             optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP),             optional, intent(in)    :: beta   !< coeff of vb [1]

  !$omp master
  if (size(v,4) /= this%ne + this%ng) then
    call Error( 'ToGhost_Merge_A'                 &
              , 'size(v,4) /= this%ne + this%ng' &
              , 'CART__Element_Transfer_Buffer'  )
  end if
  !$omp end master

  call ToGhost_Merge_X(this, v, alpha, beta)

end subroutine ToGhost_Merge_A

!-------------------------------------------------------------------------------
!> Executes MPI_Waitall to complete send of master element data

subroutine ToGhost_Finish(this)
  class(ElementTransferBuffer), intent(inout) :: this  !< buffer

  type(MPI_Status), allocatable :: status(:)
  integer :: n

  !$omp master

  n = size(this%master%request)

  if (n > 0) then
    allocate(status(n))
    call MPI_Waitall(n, this%master%request, status)
  end if

  !$omp end master
  !$omp barrier

end subroutine ToGhost_Finish

!===============================================================================

!-------------------------------------------------------------------------------
!> Send ghost data to master and receive data from own ghosts -- eXplicit

subroutine ToMaster_Transfer_X(this, mesh, v, tag)

  ! arguments ..................................................................

  class(ElementTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition), intent(in) :: mesh  !< mesh partition
  real(RNP),           intent(in) :: v     !< mesh variable
  integer,             intent(in) :: tag   !< message tag

  dimension :: v(this%np(1), this%np(2), this%np(3), this%ne+this%ng, this%nc)

  ! internal data ..............................................................

  integer :: dest, source
  integer :: i, l, m

  associate(master => this%master, ghost => this%ghost, comm => mesh%comm)

    ! copy ghost data into buffer ..............................................

    call CopyToBuffer( nn   = size(ghost%node)   &
                     , nm   = size(v) / this%nc  &
                     , nc   = this%nc            &
                     , node = ghost%node         &
                     , v    = v                  &
                     , vb   = ghost%buf          )

    ! send ghost data ..........................................................

    !$omp master
    associate(buf => ghost%buf, request => ghost%request)
      do i = 1, size(mesh%link)
        dest = mesh % link(i) % part
        m = ghost % start(i)
        l = ghost % len(i)
        if (l < 1) cycle
        call MPI_Isend(buf(m:), l, MPI_REAL_RNP, dest, tag, comm, request(i))
      end do
    end associate
    !$omp end master

    ! receive ghost data .......................................................

    !$omp master
    associate(buf => master%buf, request => master%request)
      do i = 1, size(mesh%link)
        source = mesh % link(i) % part
        m = master % start(i)
        l = master % len(i)
        if (l < 1) cycle
        call MPI_Irecv(buf(m:), l, MPI_REAL_RNP, source, tag, comm, request(i))
      end do
    end associate
    !$omp end master

  end associate

contains

  subroutine CopyToBuffer(nn, nm, nc, node, v, vb)
    integer,   intent(in)  :: nn
    integer,   intent(in)  :: nm
    integer,   intent(in)  :: nc
    integer,   intent(in)  :: node(nn)
    real(RNP), intent(in)  :: v(nm,nc)
    real(RNP), intent(out) :: vb(nn,nc)

    integer :: i, j

    !$omp do collapse(2) private(i,j)
    !$acc parallel loop collapse(2) present(v) copyin(node) copyout(vb)
    do j = 1, nc
    do i = 1, nn
      vb(i,j) = v(node(i), j)
    end do
    end do

  end subroutine CopyToBuffer

end subroutine ToMaster_Transfer_X

!-------------------------------------------------------------------------------
!> Send ghost data to master and receive data from own ghosts -- single
!>
!> As this operation extracts data from the ghost elements, the mesh variable
!> must be dimensioned as `v(np(1), np(2), np(3), ne+ng)`, where
!>
!>    *  `np` is the number of points per direction, as in `this%np`
!>    *  `ne` is the number of master elements
!>    *  `ng` is the number of ghost elements

subroutine ToMaster_Transfer_S(this, mesh, v, tag)

  ! arguments ..................................................................

  class(ElementTransferBuffer), asynchronous, intent(inout) :: this !< buffer
  type(MeshPartition),           intent(in) :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
  integer,                       intent(in) :: tag   !< message tag

  !$omp master
  if (size(v,4) /= this%ne + this%ng) then
    call Error( 'ToMaster_Transfer_S'            &
              , 'size(v,4) /= this%ne + this%ng' &
              , 'CART__Element_Transfer_Buffer'  )
  end if
  !$omp end master

  call ToMaster_Transfer_X(this, mesh, v, tag)

end subroutine ToMaster_Transfer_S

!-------------------------------------------------------------------------------
!> Send ghost data to master and receive data from own ghosts -- array
!>
!> As this operation extracts data from the ghost elements, the mesh variable
!> must be dimensioned as `v(np(1), np(2), np(3), ne+ng, nc)`, where
!>
!>    *  `np` is the number of points per direction, as in `this%np`
!>    *  `ne` is the number of master elements
!>    *  `ng` is the number of ghost elements
!>    *  `nc` is the number of components, as in `this%nc`

subroutine ToMaster_Transfer_A(this, mesh, v, tag)

  ! arguments ..................................................................

  class(ElementTransferBuffer), asynchronous,intent(inout) :: this !< buffer
  type(MeshPartition),             intent(in) :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
  integer,                         intent(in) :: tag   !< message tag

  !$omp master
  if (size(v,4) /= this%ne + this%ng) then
    call Error( 'ToMaster_Transfer_A'            &
              , 'size(v,4) /= this%ne + this%ng' &
              , 'CART__Element_Transfer_Buffer'  )
  end if
  !$omp end master

  call ToMaster_Transfer_X(this, mesh, v, tag)

end subroutine ToMaster_Transfer_A

!-------------------------------------------------------------------------------
!> Complete receive and merge buffer into into master data -- eXplicit
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>      v  =  alpha * v  +  beta * vb

subroutine ToMaster_Merge_X(this, ne_t, v, alpha, beta)

  ! arguments ..................................................................

  class(ElementTransferBuffer), intent(inout) :: this   !< buffer
  integer,                      intent(in)    :: ne_t   !< total num elements
  real(RNP),                    intent(inout) :: v      !< mesh variable
  real(RNP),          optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP),          optional, intent(in)    :: beta   !< coeff of vb [1]

  dimension :: v(this%np(1), this%np(2), this%np(3), ne_t, this%nc)

  type(MPI_Status), allocatable :: status(:)

  real(RNP) :: a, b
  integer   :: i, n

  ! initialization .............................................................

  if (size(this % master % buf) == 0) return

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

  call ScaleArray(v, a, multi=.true.)

  ! wait for receive to complete ...............................................

  !$omp master
  n = size(this % master % request)
  allocate(status(n))
  call MPI_Waitall(n, this%master%request, status)
  !$omp end master
  !$omp barrier

  ! merge buffer ...............................................................

  ! to avoid race conditions, i.e. several ghost merging to one master, the
  ! former are processed one partition (message) after the other

  associate( start => this % master % start, &
             len   => this % master % len,   &
             node  => this % master % node,  &
             vb    => this % master % buf    )

    !$acc data copyin(node,vb)
    do i = 1, size(start)
      if (len(i) < 1) cycle
      call MergeBuffer( nn   = len(i)             &
                      , nm   = size(v) / this%nc  &
                      , nc   = this%nc            &
                      , node = node(start(i):)    &
                      , b    = b                  &
                      , vb   = vb(start(i):)      &
                      , v    = v                  )
    end do
    !$acc end data

  end associate

contains

  subroutine MergeBuffer(nn, nm, nc, node, b, vb, v)
    integer,   intent(in)    :: nn
    integer,   intent(in)    :: nm
    integer,   intent(in)    :: nc
    integer,   intent(in)    :: node(nn)
    real(RNP), intent(in)    :: b
    real(RNP), intent(in)    :: vb(nn,nc)
    real(RNP), intent(inout) :: v(nm,nc)

    integer :: i, j

    !$omp do collapse(2) private(i,j)
    !$acc parallel loop collapse(2) present(v) copyin(node,vb)
    do j = 1, nc
    do i = 1, nn
      v(node(i), j) = v(node(i), j)  +  b * vb(i,j)
    end do
    end do

  end subroutine MergeBuffer

end subroutine ToMaster_Merge_X

!-------------------------------------------------------------------------------
!> Complete receive and merge buffer into ghost data -- single variable
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>      v  =  alpha * v  +  beta * vb

subroutine ToMaster_Merge_S(this, v, alpha, beta)
  class(ElementTransferBuffer),  intent(inout) :: this   !< buffer
  real(RNP), dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
  real(RNP),           optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP),           optional, intent(in)    :: beta   !< coeff of vb [1]

  call ToMaster_Merge_X(this, size(v,4), v, alpha, beta)

end subroutine ToMaster_Merge_S

!-------------------------------------------------------------------------------
!> Complete receive and merge buffer into ghost data -- array variable
!>
!> Denoting the buffer with `vb`, the following operation will be executed:
!>
!>      v  =  alpha * v  +  beta * vb

subroutine ToMaster_Merge_A(this, v, alpha, beta)
  class(ElementTransferBuffer),    intent(inout) :: this   !< buffer
  real(RNP), dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
  real(RNP),             optional, intent(in)    :: alpha  !< coeff of v  [1]
  real(RNP),             optional, intent(in)    :: beta   !< coeff of vb [1]

  call ToMaster_Merge_X(this, size(v,4), v, alpha, beta)

end subroutine ToMaster_Merge_A

!-------------------------------------------------------------------------------
!> Executes MPI_Waitall to complete send of ghost element data

subroutine ToMaster_Finish(this)
  class(ElementTransferBuffer), intent(inout) :: this  !< buffer

  type(MPI_Status), allocatable :: status(:)
  integer :: n

  !$omp master

  n = size(this%ghost%request)

  if (n > 0) then
    allocate(status(n))
    call MPI_Waitall(n, this%ghost%request, status)
  end if

  !$omp end master
  !$omp barrier

end subroutine ToMaster_Finish

!===============================================================================

!-------------------------------------------------------------------------------
!> Delete an element transfer buffer

subroutine Delete_TransferBuffer(this)
  type(ElementTransferBuffer), intent(inout) :: this  !< buffer

  this % np = 0
  this % ne = 0
  this % ng = 0
  this % nc = 0
  this % nl = 0

end subroutine Delete_TransferBuffer

!===============================================================================

end module CART__Element_Transfer_Buffer
