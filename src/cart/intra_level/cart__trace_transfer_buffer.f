!> summary:  Type for transferring linked traces (double-valued face data)
!> author:   Joerg Stiller
!> date:     2018/03/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Type for transferring linked traces (double-valued face data)
!===============================================================================

module CART__Trace_Transfer_Buffer

  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  use XMPI
  use CART__Mesh_Partition

  implicit none
  private

  public :: TraceTransferBuffer

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked face data and metadata

  type TraceTransferData
    integer,           allocatable :: side(:)    !< side to transfer (1 or 2)
    real(RNP),         allocatable :: buf(:)     !< message buffer
    type(MPI_Request), allocatable :: request(:) !< requests
  contains
    final :: Delete_TraceTransferData
  end type TraceTransferData

  !-----------------------------------------------------------------------------
  !> Type for transferring linked face data
  !>
  !> The type supports two flavors of trace variables
  !>
  !>    *  single variable of the shape `v(np,np,2,nf)`
  !>    *  variable arrays of the shape `v(np,np,2,nf,nc)`
  !>
  !>  where `v` is of type `real(RNP)` and
  !>
  !>    *  `nf` is the number of local mesh faces
  !>    *  `nc` is the number of components (variables)
  !>
  !>  Note that `v` is double-valued: `v(:,:,k,f)` is the trace of element
  !>  `face(f)%element(k)` on side `k=1` or `2` of face `f`.
  !>
  !>  Use with single thread (no OpenMP):
  !>
  !>        type(TraceTransferBuffer) :: trace_buf
  !>        ...
  !>        ! create and fill buffer, start transfer
  !>        call trace_buf % Transfer(mesh, v, tag)
  !>        ...
  !>        ! possibly perform some computations to hide communication costs
  !>        ...
  !>        ! merge received data and wait transfer to finish
  !>        call trace_buf % Merge(mesh, v, alpha, beta)
  !>        call trace_buf % Finish()
  !>
  !> Use with OpenMP:
  !>
  !>    *  the buffer must be declared `save` to be shared among the threads
  !>    *  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>
  !>          type(TraceTransferBuffer), allocatable, save :: trace_buf
  !>          ...
  !>          !$omp single
  !>          allocate(trace_buf)
  !>          !$omp end single
  !>          ...
  !>          call trace_buf % Transfer(mesh, v, tag)
  !>          ...
  !>          call trace_buf % Merge(mesh, v, alpha, beta)
  !>          call trace_buf % Finish()
  !>          ...
  !>          !$omp barrier
  !>          !$omp single
  !>          deallocate(trace_buf)
  !>          !$omp end single

  type TraceTransferBuffer
!    private

    integer :: np = 0                   !< number of points per face
    integer :: nc = 0                   !< number of components (variables)
    integer, allocatable    :: start(:) !< message start addresses
    integer, allocatable    :: len(:)   !< message lengths
    type(TraceTransferData) :: send     !< send buffer and metadata
    type(TraceTransferData) :: recv     !< receive buffer and metadata

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

    procedure :: Finish => FinishTransfer

    final :: Delete_TransferBuffer

  end type TraceTransferBuffer

contains

!===============================================================================
! TraceTransferData type bound procedures

!-------------------------------------------------------------------------------
!> Finalize trace transfer data

subroutine Delete_TraceTransferData(this)
  type(TraceTransferData), intent(inout) :: this  !< buffer

  if (allocated(this % side   )) deallocate(this % side   )
  if (allocated(this % buf    )) deallocate(this % buf    )
  if (allocated(this % request)) deallocate(this % request)

end subroutine Delete_TraceTransferData

!===============================================================================
! TraceTransferBuffer: type bound procedures

!-------------------------------------------------------------------------------
!> Create a new trace transfer buffer from mesh and given variable dimensions

subroutine New_TransferBuffer_X(this, mesh, np, nc)
  class(TraceTransferBuffer), intent(inout) :: this  !< buffer
  type(MeshPartition),        intent(in)    :: mesh  !< mesh partition
  integer,                    intent(in)    :: np    !< num points/face
  integer,                    intent(in)    :: nc    !< num components

  integer :: lb, nl, nf
  integer :: i, j

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

  ! sides .....................................................................

  nf = sum( mesh%link%nf )
  allocate(this % send % side( nf ))
  allocate(this % recv % side( nf ))

  ! sides to be sent
  associate(side => this%send%side)

    ! set side to linked face IDs
    j = 1
    do i = 1, nl
      nf = mesh%link(i)%nf
      if (nf < 1) cycle
      side(j:j+nf-1) = mesh%link(i)%face
      j = j + nf
    end do

    ! set side to position of local mesh face
    where( mesh % face(side) % element(1) >  0  .and.  &
           mesh % face(side) % element(1) <= mesh % ne )
      ! element 1 is local
      side = 1
    elsewhere
      side = 2
    end where

  end associate

  ! sides to be received
  this%recv%side = 3 - this%send%side

  ! requests ..................................................................

  allocate( this % send % request(nl), source = MPI_REQUEST_NULL )
  allocate( this % recv % request(nl), source = MPI_REQUEST_NULL )

  !$omp end single

end subroutine New_TransferBuffer_X

!-------------------------------------------------------------------------------
!> Create a new trace transfer buffer for a single variable

subroutine New_TransferBuffer_S(this, mesh, v)
  class(TraceTransferBuffer), intent(inout) :: this        !< buffer
  type(MeshPartition),        intent(in)    :: mesh        !< mesh partition
  real(RNP),                  intent(in)    :: v(:,:,:,:)  !< trace variable

  integer :: np

  np = size(v,1) * size(v,2)

  call New_TransferBuffer_X(this, mesh, np, nc=1)

end subroutine New_TransferBuffer_S

!-------------------------------------------------------------------------------
!> Create a new trace transfer buffer for an array of variables

subroutine New_TransferBuffer_A(this, mesh, v)
  class(TraceTransferBuffer), intent(inout) :: this          !< buffer
  type(MeshPartition),        intent(in)    :: mesh          !< mesh partition
  real(RNP),                  intent(in)    :: v(:,:,:,:,:)  !< trace variable

  integer :: np, nc

  np = size(v,1) * size(v,2)
  nc = size(v,5)

  call New_TransferBuffer_X(this, mesh, np, nc)

end subroutine New_TransferBuffer_A

!-------------------------------------------------------------------------------
!> Extract and transfer buffer -- eXplicit shape version

subroutine Transfer_X(this, mesh, v, tag)
  class(TraceTransferBuffer), intent(inout) :: this        !< buffer
  type(MeshPartition), intent(in) :: mesh                  !< mesh partition
  real(RNP), intent(in) :: v(this%np, 2, mesh%nf, this%nc) !< trace variable
  integer,   intent(in) :: tag                             !< message tag

  integer :: i, l, m, part
  integer :: b1, s1

  associate(send => this%send, recv => this%recv, comm => mesh%comm)

    if (size(send%buf) < 1) return

    ! extract send buffer ......................................................

    s1 = 1  ! first entry in send%side
    b1 = 1  ! first entry in send%buf

    do l = 1, size(mesh%link)

      if (mesh%link(l)%part /= mesh%part) then

        call CopyToBuffer( np   = this%np            &
                         , nc   = this%nc            &
                         , nf   = mesh%link(l)%nf    &
                         , nm   = mesh%nf            &
                         , face = mesh%link(l)%face  &
                         , side = send%side(s1:)     &
                         , v    = v                  &
                         , vb   = send%buf(b1:)      )

      else

        ! local link: copy coupled face data into send buffer,
        ! note that for the coupled face the send side equals
        ! the recv side of the linked face

        call CopyToBuffer( np   = this%np                    &
                         , nc   = this%nc                    &
                         , nf   = mesh%link(l)%nf            &
                         , nm   = mesh%nf                    &
                         , face = mesh%link(l)%coupled_face  &
                         , side = recv%side(s1:)             &
                         , v    = v                          &
                         , vb   = send%buf(b1:)              )
      end if


      s1 = s1 + mesh%link(l)%nf
      b1 = b1 + mesh%link(l)%nf * this%np * this%nc

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
                        send%request(i)                                 )

        call MPI_Irecv( recv%buf(m:), l, MPI_REAL_RNP, part, tag, comm, &
                        recv%request(i)                                 )

      else
        recv%buf(m:m+l-1) = send%buf(m:m+l-1)
      end if

    end do

    !$omp end master

  end associate

contains

  subroutine CopyToBuffer(np, nc, nf, nm, face, side, v, vb)
    integer,   intent(in)    :: np             !< number of points/face
    integer,   intent(in)    :: nc             !< number of components
    integer,   intent(in)    :: nf             !< number of linked faces
    integer,   intent(in)    :: nm             !< number of mesh faces
    integer,   intent(in)    :: face(nf)       !< list of linked faces
    integer,   intent(in)    :: side(nf)       !< linked face sides
    real(RNP), intent(in)    :: v(np,2,nm,nc)  !< trace variable
    real(RNP), intent(inout) :: vb(np*nf*nc)   !< buffer

    integer :: i, j, k

    !$omp do collapse(2)
    !$acc parallel loop collapse(3) present(v) copyin(face,side) copyout(vb)
    do k = 1, nc
    do j = 1, nf
    do i = 1, np
      vb(i + np*(j-1 + nf*(k-1))) = v(i, side(j), face(j), k)
    end do
    end do
    end do

  end subroutine CopyToBuffer

end subroutine Transfer_X

!-------------------------------------------------------------------------------
!> Extract and transfer buffer -- single variable
!>
!> The variable must be dimensioned `v(n1,n2,2,mesh%nf)` with `n1*n2 = this%np`

subroutine Transfer_S(this, mesh, v, tag)
  class(TraceTransferBuffer), intent(inout) :: this       !< buffer
  type(MeshPartition),        intent(in)    :: mesh       !< mesh partition
  real(RNP),                  intent(in)    :: v(:,:,:,:) !< trace variable
  integer,                    intent(in)    :: tag        !< message tag

  if (this%np /= size(v,1)*size(v,2) .or. this%nc /= 1) then
    call this % New(mesh, v)
  end if

  call Transfer_X(this, mesh, v, tag)

end subroutine Transfer_S

!-------------------------------------------------------------------------------
!> Extract and transfer buffer -- array of variables
!>
!> The variable must be dimensioned `v(n1,n2,2,mesh%nf,this%nc)`
!> with `n1*n2 = this%np`

subroutine Transfer_A(this, mesh, v, tag)
  class(TraceTransferBuffer), intent(inout) :: this         !< buffer
  type(MeshPartition),        intent(in)    :: mesh         !< mesh partition
  real(RNP),                  intent(in)    :: v(:,:,:,:,:) !< trace variable
  integer,                    intent(in)    :: tag          !< message tag

  if (this%np /= size(v,1)*size(v,2) .or. this%nc /= size(v,5)) then
    call this % New(mesh, v)
  end if

  call Transfer_X(this, mesh, v, tag)

end subroutine Transfer_A

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge buffer into trace variable
!>
!> Denoting the buffer with vb, the following operation will be executed:
!>
!>    `v  =  alpha * v  +  beta * vb`

subroutine Merge_X(this, mesh, v, alpha, beta)
  class(TraceTransferBuffer), intent(inout) :: this           !< buffer
  type(MeshPartition), intent(in) :: mesh                     !< mesh partition
  real(RNP), intent(inout) :: v(this%np, 2, mesh%nf, this%nc) !< trace variable
  real(RNP), optional, intent(in) :: alpha                    !< coeff of v  [1]
  real(RNP), optional, intent(in) :: beta                     !< coeff of vb [1]

  type(MPI_Status), allocatable :: status(:)

  integer :: l, nr
  integer :: b1, s1
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

  s1 = 1  ! first entry in recv%side
  b1 = 1  ! first entry in recv%buf

  do l = 1, size(mesh%link)

    call MergeBuffer(  np   = this%np              &
                    ,  nc   = this%nc              &
                    ,  nf   = mesh%link(l)%nf      &
                    ,  nm   = mesh%nf              &
                    ,  face = mesh%link(l)%face    &
                    ,  side = this%recv%side(s1:)  &
                    ,  a    = a                    &
                    ,  b    = b                    &
                    ,  v    = v                    &
                    ,  vb   = this%recv%buf(b1:)   )

    s1 = s1 + mesh%link(l)%nf
    b1 = b1 + mesh%link(l)%nf * this%np * this%nc

  end do

contains

  subroutine MergeBuffer(np, nc, nf, nm, face, side, a, b, v, vb)
    integer,   intent(in)    :: np             !< number of points/face
    integer,   intent(in)    :: nc             !< number of components
    integer,   intent(in)    :: nf             !< number of linked faces
    integer,   intent(in)    :: nm             !< number of mesh faces
    integer,   intent(in)    :: face(nf)       !< list of linked faces
    integer,   intent(in)    :: side(nf)       !< linked face sides
    real(RNP), intent(in)    :: a              !< coefficient of v
    real(RNP), intent(in)    :: b              !< coefficient of vb
    real(RNP), intent(inout) :: v(np,2,nm,nc)  !< face variable
    real(RNP), intent(in)    :: vb(np*nf*nc)   !< buffer

    integer :: i, j, k

    !$omp do collapse(2)
    !$acc parallel loop collapse(3) present(v) copyin(face,side,vb)
    do k = 1, nc
    do j = 1, nf
    do i = 1, np
      v(i, side(j), face(j), k)  =  a * v(i, side(j), face(j), k)  &
                                 +  b * vb(i + np*(j-1 + nf*(k-1)))
    end do
    end do
    end do

  end subroutine MergeBuffer

end subroutine Merge_X

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge real buffer -- single variable
!>
!> The variable must be dimensioned `v(n1,n2,2,mesh%nf)` with `n1*n2 = this%np`

subroutine Merge_S(this, mesh, v, alpha, beta)
  class(TraceTransferBuffer), intent(inout) :: this        !< buffer
  type(MeshPartition),        intent(in)    :: mesh        !< mesh partition
  real(RNP),                  intent(inout) :: v(:,:,:,:)  !< face variable
  real(RNP),        optional, intent(in)    :: alpha       !< coeff of v  [1]
  real(RNP),        optional, intent(in)    :: beta        !< coeff of vb [1]

  call Merge_X(this, mesh, v, alpha, beta)

end subroutine Merge_S

!-------------------------------------------------------------------------------
!> Wait receive to complete and merge real buffer -- array of variables
!>
!> The variable must be dimensioned `v(n1,n2,2,mesh%nf,this%nc)`
!> with `n1*n2 = this%np`

subroutine Merge_A(this, mesh, v, alpha, beta)
  class(TraceTransferBuffer), intent(inout) :: this          !< buffer
  type(MeshPartition),        intent(in)    :: mesh          !< mesh partition
  real(RNP),                  intent(inout) :: v(:,:,:,:,:)  !< face variable
  real(RNP),        optional, intent(in)    :: alpha         !< coeff of v  [1]
  real(RNP),        optional, intent(in)    :: beta          !< coeff of vb [1]

  call Merge_X(this, mesh, v, alpha, beta)

end subroutine Merge_A

!-------------------------------------------------------------------------------
!> Executes MPI_Waitall to complete all transfers associated with the buffer

subroutine FinishTransfer(this)
  class(TraceTransferBuffer), intent(inout) :: this  !< buffer

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
!> Finalization of a trace transfer buffer

subroutine Delete_TransferBuffer(this)
  type(TraceTransferBuffer), intent(inout) :: this  !< buffer

  this % np = 0
  this % nc = 0

  if (allocated(this % start)) deallocate(this % start)
  if (allocated(this % len  )) deallocate(this % len  )

end subroutine Delete_TransferBuffer

!===============================================================================

end module CART__Trace_Transfer_Buffer
