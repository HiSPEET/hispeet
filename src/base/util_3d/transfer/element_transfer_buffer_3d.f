!> summary:  Type and methods for transferring linked element data
!> author:   Joerg Stiller
!> date:     2017/09/05, revised 2021/01/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!> @todo
!>   - check wether OpenMP parallelization makes sense
!===============================================================================

module Element_Transfer_Buffer_3d

  use Kind_Parameters,   only: RNP
! use Array_Assignments, only: ScaleArray
  use Execution_Control, only: Error
  use XMPI
  use Mesh_3d__Link
  use Mesh_3d__Partition

  implicit none
  private

  public :: ElementTransferBuffer3d

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked element data and metadata

  type ElementTransferData
    integer,           allocatable :: start(:)   !< message start addresses
    integer,           allocatable :: len(:)     !< message lengths
    integer,           allocatable :: node(:)    !< mesh point IDs
    integer,           allocatable :: ibuf(:)    !< integer message buffer
    real(RNP),         allocatable :: rbuf(:)    !< real message buffer
    type(MPI_Request), allocatable :: request(:) !< MPI requests
  contains
    final :: Delete_ElementTransferData
  end type ElementTransferData

  !-----------------------------------------------------------------------------
  !> Type for transferring data from master to ghost elements
  !>
  !> For realizing anisotropic, element-overlapping operations, the number
  !> of element "points", `np`, can vary in each coordinate direction.
  !> Moreover, the implementation supports full and partial transfers.
  !> In the latter case, only `nl` point layers are transferred. The number
  !> of layers can be chosen independently for each coordinate direction.
  !> If `nl` is not specified or identical to zero, a full copy is performed.
  !>
  !> The type supports two flavors of mesh variables
  !>
  !>    -  single variable of the shape `v(np(1),np(2),np(3),ne+ng)`
  !>    -  variable arrays of the shape `v(np(1),np(2),np(3),ne+ng,nc)`
  !>
  !> where `v` is either of type `real(RNP)` or `integer` and
  !>
  !>    -  `np(1:3)` is the number of element or region points per direction
  !>    -  `ne` is the number of local elements, i.e. `mesh%ne`
  !>    -  `ng` is the number of ghost elements, i.e. `mesh%ng`
  !>    -  `nc` is the number of components (variables)
  !>
  !> It is recommended to pass the variable arguments always in the shape
  !> given above, though the ghost entries are not used and can be omitted
  !> in certain cases.
  !>
  !> Use with single thread (no OpenMP):
  !>
  !>         type(ElementTransferBuffer3d) :: buf
  !>         ...
  !>         ! create and fill buffer, start transfer
  !>         buf = ElementTransferBuffer3d(mesh, v, nl)
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         ! possibly perform some computations to hide communication costs
  !>         ...
  !>         ! merge received master data into ghosts and finish transfer
  !>         call buf % Merge(v, alpha, beta)
  !>         call buf % Finish()
  !>
  !> Use with OpenMP:
  !>
  !>    -  the buffer must be declared `save` to be shared among the threads
  !>    -  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>
  !>         type(ElementTransferBuffer3d), allocatable, save :: buf
  !>         ...
  !>         !$omp single
  !>         buf = ElementTransferBuffer3d(mesh, v, nl)
  !>         !$omp end single
  !>         ...
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         call buf % Merge(v, alpha, beta)
  !>         call buf % Finish()
  !>         ...
  !>         !$omp barrier
  !>         !$omp single
  !>         deallocate(buf)
  !>         !$omp end single

  type ElementTransferBuffer3d

    integer :: np(3) = 0  !< points per direction
    integer :: ne    = 0  !< number of local (master) elements
    integer :: ng    = 0  !< number of ghost elements
    integer :: nc    = 0  !< number of components (variables)
    integer :: nl(3) = 0  !< number of layers, all(nl < 1) implies full transfer

    type(ElementTransferData) :: master  !< master element data
    type(ElementTransferData) :: ghost   !< ghost element data

  contains

    ! public procedures ........................................................

    generic,   public  :: Init     => Init_IS, Init_IA, &
                                      Init_RS, Init_RA

    generic,   public  :: Transfer => Transfer_IS, Transfer_IA, &
                                      Transfer_RS, Transfer_RA

    generic,   public  :: Merge    => Merge_IS, Merge_IA, &
                                      Merge_RS, Merge_RA
    procedure, public  :: Finish

    ! private procedures .......................................................

    procedure, private :: Init_IS, Init_IA
    procedure, private :: Init_RS, Init_RA

    procedure, private :: Transfer_IS, Transfer_IA
    procedure, private :: Transfer_RS, Transfer_RA

    procedure, private :: Merge_IS, Merge_IA
    procedure, private :: Merge_RS, Merge_RA

    ! finalization .............................................................

    final :: Delete_TransferBuffer

  end type ElementTransferBuffer3d

  !-----------------------------------------------------------------------------
  ! Constructor

  interface ElementTransferBuffer3d
    module procedure New_TransferBuffer_IS
    module procedure New_TransferBuffer_IA
    module procedure New_TransferBuffer_RS
    module procedure New_TransferBuffer_RA
  end interface

contains

  !=============================================================================
  ! ElementTransferData: type bound procedures

  !-----------------------------------------------------------------------------
  !> Delete element transfer data

  subroutine Delete_ElementTransferData(this)
    type(ElementTransferData), intent(inout) :: this  !< buffer

    if (allocated(this % start  )) deallocate(this % start  )
    if (allocated(this % len    )) deallocate(this % len    )
    if (allocated(this % node   )) deallocate(this % node   )
    if (allocated(this % ibuf   )) deallocate(this % ibuf   )
    if (allocated(this % rbuf   )) deallocate(this % rbuf   )
    if (allocated(this % request)) deallocate(this % request)

  end subroutine Delete_ElementTransferData

  !=============================================================================
  ! ElementTransferBuffer3d: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single integer variable

  function New_TransferBuffer_IS(mesh, v, nl) result(this)
    type(Mesh3d_Partition),      intent(in) :: mesh  !< mesh partition
    integer, dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,           optional, intent(in) :: nl(3) !< number of layers
    type(ElementTransferBuffer3d) :: this

    call Init_IS(this, mesh, v, nl)

  end function New_TransferBuffer_IS

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of integer variables

  function New_TransferBuffer_IA(mesh, v, nl) result(this)
    type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
    integer, dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,             optional, intent(in) :: nl(3) !< number of layers
    type(ElementTransferBuffer3d) :: this

    call Init_IA(this, mesh, v, nl)

  end function New_TransferBuffer_IA

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single real variable

  function New_TransferBuffer_RS(mesh, v, nl) result(this)
    type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,             optional, intent(in) :: nl(3) !< number of layers
    type(ElementTransferBuffer3d) :: this

    call Init_RS(this, mesh, v, nl)

  end function New_TransferBuffer_RS

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of real variables

  function New_TransferBuffer_RA(mesh, v, nl) result(this)
    type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,               optional, intent(in) :: nl(3) !< number of layers
    type(ElementTransferBuffer3d) :: this

    call Init_RA(this, mesh, v, nl)

  end function New_TransferBuffer_RA

  !=============================================================================
  ! ElementTransferBuffer3d: initialization procedures

  !-----------------------------------------------------------------------------
  !> Create a new element transfer buffer for a single integer variable

  subroutine Init_IS(this, mesh, v, nl)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    integer,    dimension(:,:,:,:), intent(in)    :: v     !< mesh variable
    integer,              optional, intent(in)    :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_X(this, mesh, np, nl)

    this % nc = 1

    allocate(this % master % ibuf( size(this % master % node) ))
    allocate(this % ghost  % ibuf( size(this % ghost  % node) ))

  end subroutine Init_IS

  !-----------------------------------------------------------------------------
  !> Create a new element transfer buffer for an array of integer variables

  subroutine Init_IA(this, mesh, v, nl)
    class(ElementTransferBuffer3d),  intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),          intent(in)    :: mesh  !< mesh partition
    integer,   dimension(:,:,:,:,:), intent(in)    :: v     !< mesh variables
    integer,               optional, intent(in)    :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_X(this, mesh, np, nl)

    this % nc = nc

    allocate(this % master % ibuf( size(this % master % node) * nc ))
    allocate(this % ghost  % ibuf( size(this % ghost  % node) * nc ))

  end subroutine Init_IA

  !-----------------------------------------------------------------------------
  !> Create a new element transfer buffer for a single real variable

  subroutine Init_RS(this, mesh, v, nl)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    real(RNP),  dimension(:,:,:,:), intent(in)    :: v     !< mesh variable
    integer,              optional, intent(in)    :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_X(this, mesh, np, nl)

    this % nc = 1

    allocate(this % master % rbuf( size(this % master % node) ))
    allocate(this % ghost  % rbuf( size(this % ghost  % node) ))

  end subroutine Init_RS

  !-----------------------------------------------------------------------------
  !> Create a new element transfer buffer for an array of real variables

  subroutine Init_RA(this, mesh, v, nl)
    class(ElementTransferBuffer3d),  intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),          intent(in)    :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:,:), intent(in)    :: v     !< mesh variables
    integer,               optional, intent(in)    :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_X(this, mesh, np, nl)

    this % nc = nc

    allocate(this % master % rbuf( size(this % master % node) * nc ))
    allocate(this % ghost  % rbuf( size(this % ghost  % node) * nc ))

  end subroutine Init_RA

  !-----------------------------------------------------------------------------
  !> Common initialization of element transfer buffers

  subroutine Init_X(this, mesh, np, nl)
    class(ElementTransferBuffer3d), intent(inout) :: this !< buffer
    type(Mesh3d_Partition), intent(in) :: mesh  !< mesh partition
    integer,                intent(in) :: np(3) !< points per direction
    integer,      optional, intent(in) :: nl(3) !< number of layers

    logical, allocatable :: mask(:,:,:)
    integer, allocatable :: node(:)

    integer :: lm, lg, np1, np12, np123
    integer :: e, i, j, k, l, m, n
    logical :: partial

    !$omp single

    ! setup ....................................................................

    this % np = np
    this % ne = mesh % n_elem
    this % ng = mesh % n_ghost
    if (present(nl)) then
      this % nl = nl
    else
      this % nl = 0
    end if
    partial = any(this%nl > 0)

    lm = sum(mesh % link % n_master)
    lg = sum(mesh % link % n_ghost )

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
            do k = 1, np(1)
            do j = 1, np(2)
            do i = 1, np(3)
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

      allocate( this % master % node(n), source = node(1:n) )
      allocate( this % master % request(size(mesh%link)), &
                source = MPI_REQUEST_NULL )

    end associate

    ! ghost element data .......................................................

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
            do k = 1, np(1)
            do j = 1, np(2)
            do i = 1, np(3)
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

      allocate( this % ghost % node(n), source = node(1:n) )

      allocate( this % ghost % request(size(mesh%link)), &
                source = MPI_REQUEST_NULL )

    end associate

    !$omp end single

  end subroutine Init_X

  !-------------------------------------------------------------------------------
  !> Generates a mask of points subjected to transfer operations

  pure subroutine GeneratePointMask(link, np, nl, mask)
    class(Mesh3d_ElementLink), intent(in)  :: link !< element link
    integer, intent(in)  :: np(3)       !< points per direction
    integer, intent(in)  :: nl(3)       !< point layers per direction
    logical, intent(out) :: mask(:,:,:) !< mask of linked points

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

  !=============================================================================
  ! ElementTransferBuffer3d: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer scalar

  subroutine Transfer_IS(this, mesh, v, tag)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    integer, dimension(:,:,:,:),    intent(in)    :: v     !< mesh variable
    integer,                        intent(in)    :: tag   !< message tag

    call Transfer_IX(this, mesh, size(v,4), v, tag)

  end subroutine Transfer_IS

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer array

  subroutine Transfer_IA(this, mesh, v, tag)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    integer, dimension(:,:,:,:,:),  intent(in)    :: v     !< mesh variable
    integer,                        intent(in)    :: tag   !< message tag

    call Transfer_IX(this, mesh, size(v,4), v, tag)

  end subroutine Transfer_IA

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer eXplicit

  subroutine Transfer_IX(this, mesh, ne_t, v, tag)

    ! arguments ................................................................

    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    integer,                        intent(in)    :: ne_t  !< total num elements
    integer,                        intent(in)    :: v     !< mesh variable
    integer,                        intent(in)    :: tag   !< message tag

    dimension :: v(this%np(1), this%np(2), this%np(3), ne_t, this%nc)

    ! internal data ..............................................................

    integer :: dest, source
    integer :: i, l, m

    ! sanity check ..............................................................

    !$omp master
    if (ne_t < this%ne) then
      call Error('Transfer_IX','size(v,4) < this%ne','Element_Transfer_Buffer_3d')
    end if
    !$omp end master

    associate(master => this%master, ghost => this%ghost, comm => mesh%comm)

      ! copy master data into buffer ............................................

      call CopyToBuffer( nn   = size(master%node)  &
                       , nm   = size(v) / this%nc  &
                       , nc   = this%nc            &
                       , node = master%node        &
                       , v    = v                  &
                       , vb   = master%ibuf        )

      ! send master data ......................................................

      !$omp master
      associate(buf => master%ibuf, request => master%request)
        do i = 1, size(mesh%link)
          dest = mesh % link(i) % part
          m = master % start(i)
          l = master % len(i)
          if (l < 1) cycle
          call MPI_Isend(buf(m:), l, MPI_INTEGER, dest, tag, comm, request(i))
        end do
      end associate
      !$omp end master

      ! receive master data ....................................................

      !$omp master
      associate(buf => ghost%ibuf, request => ghost%request)
        do i = 1, size(mesh%link)
          source = mesh % link(i) % part
          m = ghost % start(i)
          l = ghost % len(i)
          if (l < 1) cycle
          call MPI_Irecv(buf(m:), l, MPI_INTEGER, source, tag, comm, request(i))
        end do
      end associate
      !$omp end master

    end associate

  contains

    subroutine CopyToBuffer(nn, nm, nc, node, v, vb)
      integer, intent(in)  :: nn
      integer, intent(in)  :: nm
      integer, intent(in)  :: nc
      integer, intent(in)  :: node(nn)
      integer, intent(in)  :: v(nm,nc)
      integer, intent(out) :: vb(nn,nc)

      integer :: i, j

      !$omp do collapse(2)
      do j = 1, nc
      do i = 1, nn
        vb(i,j) = v(node(i), j)
      end do
      end do

    end subroutine CopyToBuffer

  end subroutine Transfer_IX

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real scalar

  subroutine Transfer_RS(this, mesh, v, tag)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:),  intent(in)    :: v     !< mesh variable
    integer,                        intent(in)    :: tag   !< message tag

    call Transfer_RX(this, mesh, size(v,4), v, tag)

  end subroutine Transfer_RS

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real array

  subroutine Transfer_RA(this, mesh, v, tag)
    class(ElementTransferBuffer3d),  intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),          intent(in)    :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:,:), intent(in)    :: v     !< mesh variable
    integer,                         intent(in)    :: tag   !< message tag

    call Transfer_RX(this, mesh, size(v,4), v, tag)

  end subroutine Transfer_RA

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real eXplicit

  subroutine Transfer_RX(this, mesh, ne_t, v, tag)

    ! arguments ................................................................

    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer
    type(Mesh3d_Partition),         intent(in)    :: mesh  !< mesh partition
    integer,                        intent(in)    :: ne_t  !< total num elements
    real(RNP),                      intent(in)    :: v     !< mesh variable
    integer,                        intent(in)    :: tag   !< message tag

    dimension :: v(this%np(1), this%np(2), this%np(3), ne_t, this%nc)

    ! internal data ..............................................................

    integer :: dest, source
    integer :: i, l, m

    ! sanity check ..............................................................

    !$omp master
    if (ne_t < this%ne) then
      call Error('Transfer_RX','size(v,4) < this%ne','Element_Transfer_Buffer_3d')
    end if
    !$omp end master

    associate(master => this%master, ghost => this%ghost, comm => mesh%comm)

      ! copy master data into buffer ............................................

      call CopyToBuffer( nn   = size(master%node)  &
                       , nm   = size(v) / this%nc  &
                       , nc   = this%nc            &
                       , node = master%node        &
                       , v    = v                  &
                       , vb   = master%rbuf        )

      ! send master data ......................................................

      !$omp master
      associate(buf => master%rbuf, request => master%request)
        do i = 1, size(mesh%link)
          dest = mesh % link(i) % part
          m = master % start(i)
          l = master % len(i)
          if (l < 1) cycle
          call MPI_Isend(buf(m:), l, MPI_REAL_RNP, dest, tag, comm, request(i))
        end do
      end associate
      !$omp end master

      ! receive master data ....................................................

      !$omp master
      associate(buf => ghost%rbuf, request => ghost%request)
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

      !$omp do collapse(2)
      do j = 1, nc
      do i = 1, nn
        vb(i,j) = v(node(i), j)
      end do
      end do

    end subroutine CopyToBuffer

  end subroutine Transfer_RX

  !=============================================================================
  ! ElementTransferBuffer3d: merge procedures

  !-----------------------------------------------------------------------------
  !> Complete receive merge buffer into ghost data -- integer scalar
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

  subroutine Merge_IS(this, v, alpha, beta)
    class(ElementTransferBuffer3d), intent(inout) :: this   !< buffer
    integer,    dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
    integer,              optional, intent(in)    :: alpha  !< coeff of v  [1]
    integer,              optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IS'                       &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffe_3dr'     )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IS

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer array
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

  subroutine Merge_IA(this, v, alpha, beta)
    class(ElementTransferBuffer3d),  intent(inout) :: this   !< buffer
    integer,   dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
    integer,               optional, intent(in)    :: alpha  !< coeff of v  [1]
    integer,               optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IA'                       &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IA

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer eXplicit
  !>
  !> Denoting the buffer with `vb`, the following operation will be executed:
  !>
  !>       v  =  alpha * v  +  beta * vb

  subroutine Merge_IX(this, v, alpha, beta)

    ! arguments ................................................................

    class(ElementTransferBuffer3d), intent(inout) :: this   !< buffer
    integer,                        intent(inout) :: v      !< mesh variable
    integer,              optional, intent(in)    :: alpha  !< coeff of v  [1]
    integer,              optional, intent(in)    :: beta   !< coeff of vb [1]

    dimension :: v(this%np(1), this%np(2), this%np(3), this%ne+this%ng, this%nc)

    ! internal data ............................................................

    type(MPI_Status), allocatable :: status(:)

    integer :: a, b
    integer :: n

    ! initialization ...........................................................

    if (size(this % ghost % ibuf) == 0) return

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

    ! wait for receive to complete .............................................

    !$omp master
    n = size(this % ghost % request)
    allocate(status(n))
    call MPI_Waitall(n, this%ghost%request, status)
    !$omp end master
    !$omp barrier

    ! merge buffer .............................................................

    call MergeBuffer( nn   = size(this%ghost%node)  &
                    , nm   = size(v) / this % nc    &
                    , nc   = this % nc              &
                    , node = this % ghost % node    &
                    , a    = a                      &
                    , b    = b                      &
                    , vb   = this % ghost % ibuf    &
                    , v    = v                      )

  contains

    subroutine MergeBuffer(nn, nm, nc, node, a, b, vb, v)
      integer, intent(in)    :: nn
      integer, intent(in)    :: nm
      integer, intent(in)    :: nc
      integer, intent(in)    :: node(nn)
      integer, intent(in)    :: a, b
      integer, intent(in)    :: vb(nn,nc)
      integer, intent(inout) :: v(nm,nc)

      integer :: i, j

      if (a /= ZERO) then

        !$omp do collapse(2)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = a * v(node(i), j)  +  b * vb(i,j)
        end do
        end do

      else

        !$omp do collapse(2)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = b * vb(i,j)
        end do
        end do

      end if

    end subroutine MergeBuffer

  end subroutine Merge_IX

  !-----------------------------------------------------------------------------
  !> Complete receive merge buffer into ghost data -- real scalar
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

  subroutine Merge_RS(this, v, alpha, beta)
    class(ElementTransferBuffer3d), intent(inout) :: this   !< buffer
    real(RNP),  dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RNP),            optional, intent(in)    :: alpha  !< coeff of v  [1]
    real(RNP),            optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RS'                       &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffe_3dr'     )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RS

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real array
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

  subroutine Merge_RA(this, v, alpha, beta)
    class(ElementTransferBuffer3d),  intent(inout) :: this   !< buffer
    real(RNP), dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RNP),             optional, intent(in)    :: alpha  !< coeff of v  [1]
    real(RNP),             optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RA'                       &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RA

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real eXplicit
  !>
  !> Denoting the buffer with `vb`, the following operation will be executed:
  !>
  !>       v  =  alpha * v  +  beta * vb

  subroutine Merge_RX(this, v, alpha, beta)

    ! arguments ................................................................

    class(ElementTransferBuffer3d), intent(inout) :: this   !< buffer
    real(RNP),                      intent(inout) :: v      !< mesh variable
    real(RNP),            optional, intent(in)    :: alpha  !< coeff of v  [1]
    real(RNP),            optional, intent(in)    :: beta   !< coeff of vb [1]

    dimension :: v(this%np(1), this%np(2), this%np(3), this%ne+this%ng, this%nc)

    ! internal data ............................................................

    type(MPI_Status), allocatable :: status(:)

    real(RNP) :: a, b
    integer   :: n

    ! initialization ...........................................................

    if (size(this % ghost % rbuf) == 0) return

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

    ! wait for receive to complete .............................................

    !$omp master
    n = size(this % ghost % request)
    allocate(status(n))
    call MPI_Waitall(n, this%ghost%request, status)
    !$omp end master
    !$omp barrier

    ! merge buffer .............................................................

    call MergeBuffer( nn   = size(this%ghost%node)  &
                    , nm   = size(v) / this % nc    &
                    , nc   = this % nc              &
                    , node = this % ghost % node    &
                    , a    = a                      &
                    , b    = b                      &
                    , vb   = this % ghost % rbuf    &
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

        !$omp do collapse(2)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = a * v(node(i), j)  +  b * vb(i,j)
        end do
        end do

      else

        !$omp do collapse(2)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = b * vb(i,j)
        end do
        end do

      end if

    end subroutine MergeBuffer

  end subroutine Merge_RX

  !=============================================================================
  ! ElementTransferBuffer3d: finish procedure

  !-----------------------------------------------------------------------------
  !> Executes MPI_Waitall to complete send of master element data

  subroutine Finish(this)
    class(ElementTransferBuffer3d), intent(inout) :: this  !< buffer

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

  end subroutine Finish

  !=============================================================================
  ! ElementTransferBuffer3d: finalization

  !-----------------------------------------------------------------------------
  !> Delete element transfer buffer

  subroutine Delete_TransferBuffer(this)
    type(ElementTransferBuffer3d), intent(inout) :: this  !< buffer

    this % np = 0
    this % ne = 0
    this % ng = 0
    this % nc = 0
    this % nl = 0

  end subroutine Delete_TransferBuffer

  !=============================================================================

end module Element_Transfer_Buffer_3d
