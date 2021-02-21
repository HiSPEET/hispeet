!> summary:  Type and methods for transferring linked element data
!> author:   Joerg Stiller
!> date:     2017/09/05, revised 2021/01/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!> @todo
!>   - check wether OpenMP parallelization makes sense
!===============================================================================

module Element_Transfer_Buffer_3d

  use Kind_Parameters
  use Constants        , only: ZERO, ONE
  use Execution_Control, only: Error
  use XMPI
  use Mesh_3d__Link
  use Mesh_3d__Partition

  implicit none
  private

  public :: ElementTransferBuffer3d

  ! made public to circumvent an error with gfortran*9/10:
  public :: Init_ElementTransferBuffer3d_Shared

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked element data and metadata

  type ElementTransferData
    integer,           allocatable :: start(:)   !< message start addresses
    integer,           allocatable :: len(:)     !< message lengths
    integer,           allocatable :: node(:)    !< mesh point IDs
    class(*),          allocatable :: buf(:)     !< message buffer
    type(MPI_Request), allocatable :: request(:) !< MPI requests
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
  !> where `v` is either of type `real(RNP)`, `integer`, `integer(IXS)` or
  !> `integer(IXL)` and
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
  !>    -  the `asynchronous` attribute is required for non-blocking MPI send
  !>       and receive operations in `Transfer`
  !>
  !>         type(ElementTransferBuffer3d), asynchronous :: buf
  !>         ...
  !>         ! create and fill buffer, start transfer
  !>         buf = ElementTransferBuffer3d(mesh, v, nl)
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         ! possibly perform some computations to hide communication costs
  !>         ...
  !>         ! merge received master data into ghosts and finish transfer
  !>         call buf % Merge(v, alpha, beta)
  !>
  !> Use with OpenMP:
  !>
  !>    -  the buffer must be declared `save` to be shared among the threads
  !>    -  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>
  !>         type(ElementTransferBuffer3d), asynchronous, allocatable, save :: buf
  !>         ...
  !>         !$omp master
  !>         buf = ElementTransferBuffer3d(mesh, v, nl)
  !>         !$omp end master
  !>         !$omp barrier
  !>         ...
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         call buf % Merge(v, alpha, beta)
  !>         ...
  !>         !$omp master
  !>         deallocate(buf)
  !>         !$omp end master

  type ElementTransferBuffer3d

    integer :: np(3) = 0  !< points per direction
    integer :: ne    = 0  !< number of local (master) elements
    integer :: ng    = 0  !< number of ghost elements
    integer :: nc    = 0  !< number of components (variables)
    integer :: nl(3) = 0  !< number of layers, all(nl < 1) implies full transfer

    type(ElementTransferData) :: master  !< master element data
    type(ElementTransferData) :: ghost   !< ghost element data

  contains

    generic, public :: Transfer => Transfer_IDK_S, Transfer_IDK_A, &
                                   Transfer_IXS_S, Transfer_IXS_A, &
                                   Transfer_IXL_S, Transfer_IXL_A, &
                                   Transfer_RNP_S, Transfer_RNP_A

    generic, public :: Merge    => Merge_IDK_S, Merge_IDK_A, &
                                   Merge_IXS_S, Merge_IXS_A, &
                                   Merge_IXL_S, Merge_IXL_A, &
                                   Merge_RNP_S, Merge_RNP_A

    procedure, private :: Transfer_IDK_S, Transfer_IDK_A
    procedure, private :: Transfer_IXS_S, Transfer_IXS_A
    procedure, private :: Transfer_IXL_S, Transfer_IXL_A
    procedure, private :: Transfer_RNP_S, Transfer_RNP_A

    procedure, private :: Merge_IDK_S, Merge_IDK_A
    procedure, private :: Merge_IXS_S, Merge_IXS_A
    procedure, private :: Merge_IXL_S, Merge_IXL_A
    procedure, private :: Merge_RNP_S, Merge_RNP_A

  end type ElementTransferBuffer3d

  !=============================================================================
  ! Constructor

  interface ElementTransferBuffer3d

    !---------------------------------------------------------------------------
    !> New element transfer buffer for a single integer variable

    module function New_TransferBuffer_IDK_S(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),      intent(in) :: mesh  !< mesh partition
      integer, dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,           optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IDK_S

    !---------------------------------------------------------------------------
    !> New element transfer buffer for an array of integer variables

    module function New_TransferBuffer_IDK_A(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
      integer, dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
      integer,             optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IDK_A

    !---------------------------------------------------------------------------
    !> New element transfer buffer for a single integer(IXS) variable

    module function New_TransferBuffer_IXS_S(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),           intent(in) :: mesh  !< mesh partition
      integer(IXS), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IXS_S

    !---------------------------------------------------------------------------
    !> New element transfer buffer for an array of integer(IXS) variables

    module function New_TransferBuffer_IXS_A(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),             intent(in) :: mesh  !< mesh partition
      integer(IXS), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
      integer,                  optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IXS_A

    !---------------------------------------------------------------------------
    !> New element transfer buffer for a single integer(IXL) variable

    module function New_TransferBuffer_IXL_S(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),           intent(in) :: mesh  !< mesh partition
      integer(IXL), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IXL_S

    !---------------------------------------------------------------------------
    !> New element transfer buffer for an array of integer(IXL) variables

    module function New_TransferBuffer_IXL_A(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),             intent(in) :: mesh  !< mesh partition
      integer(IXL), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
      integer,                  optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_IXL_A

    !---------------------------------------------------------------------------
    !> New element transfer buffer for a single real variable

    module function New_TransferBuffer_RNP_S(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
      real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,             optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_RNP_S

    !---------------------------------------------------------------------------
    !> New element transfer buffer for an array of real variables

    module function New_TransferBuffer_RNP_A(mesh, v, nl) result(this)
      type(ElementTransferBuffer3d) :: this
      type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
      real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
      integer,               optional, intent(in) :: nl(3) !< number of layers
    end function New_TransferBuffer_RNP_A

  end interface

  !=============================================================================
  ! Transfer procedures

  interface

    !--------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer scalar

    module subroutine Transfer_IDK_S(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),      intent(in) :: mesh  !< mesh partition
      integer, dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                     intent(in) :: tag   !< message tag
    end subroutine Transfer_IDK_S

    !---------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer array

    module subroutine Transfer_IDK_A(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
      integer, dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                       intent(in) :: tag   !< message tag
    end subroutine Transfer_IDK_A

    !--------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer(IXS) scalar

    module subroutine Transfer_IXS_S(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),           intent(in) :: mesh  !< mesh partition
      integer(IXS), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                          intent(in) :: tag   !< message tag
    end subroutine Transfer_IXS_S

    !---------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer(IXS) array

    module subroutine Transfer_IXS_A(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),             intent(in) :: mesh  !< mesh partition
      integer(IXS), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                            intent(in) :: tag   !< message tag
    end subroutine Transfer_IXS_A

    !--------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer(IXL) scalar

    module subroutine Transfer_IXL_S(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),           intent(in) :: mesh  !< mesh partition
      integer(IXL), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                          intent(in) :: tag   !< message tag
    end subroutine Transfer_IXL_S

    !---------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- integer(IXL) array

    module subroutine Transfer_IXL_A(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),             intent(in) :: mesh  !< mesh partition
      integer(IXL), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                            intent(in) :: tag   !< message tag
    end subroutine Transfer_IXL_A

    !--------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- real scalar

    module subroutine Transfer_RNP_S(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
      real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                       intent(in) :: tag   !< message tag
    end subroutine Transfer_RNP_S

    !---------------------------------------------------------------------------
    !> Send master data to ghosts and receive ghost data -- real array

    module subroutine Transfer_RNP_A(this, mesh, v, tag)
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
      type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
      real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
      integer,                         intent(in) :: tag   !< message tag
    end subroutine Transfer_RNP_A

  end interface

  !=============================================================================
  ! Merge procedures

  interface

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer scalar

    module subroutine Merge_IDK_S(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer, dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,           optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,           optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IDK_S

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer array

    module subroutine Merge_IDK_A(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer, dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,             optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,             optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IDK_A

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer(IXS) scalar

    module subroutine Merge_IXS_S(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer(IXS), dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,                optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,                optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IXS_S

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer(IXS) array

    module subroutine Merge_IXS_A(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer(IXS), dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,                  optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,                  optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IXS_A

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer(IXL) scalar

    module subroutine Merge_IXL_S(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer(IXL), dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,                optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,                optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IXL_S

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- integer(IXL) array

    module subroutine Merge_IXL_A(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      integer(IXL), dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
      integer,                  optional, intent(in)    :: alpha !< coeff of v  [0]
      integer,                  optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_IXL_A

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- real scalar

    module subroutine Merge_RNP_S(this, v, alpha, beta)
      class(ElementTransferBuffer3d), intent(inout) :: this
      real(RNP),  dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
      real(RNP),            optional, intent(in)    :: alpha !< coeff of v  [0]
      real(RNP),            optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_RNP_S

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into ghost data -- real array

    module subroutine Merge_RNP_A(this, v, alpha, beta)
      class(ElementTransferBuffer3d) , intent(inout) :: this
      real(RNP), dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
      real(RNP),             optional, intent(in)    :: alpha !< coeff of v  [0]
      real(RNP),             optional, intent(in)    :: beta  !< coeff of vb [1]
    end subroutine Merge_RNP_A

  end interface

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Common initialization of element transfer buffers

  subroutine Init_ElementTransferBuffer3d_Shared(this, mesh, np, nl)
    class(ElementTransferBuffer3d), intent(inout) :: this !< buffer
    type(Mesh3d_Partition), intent(in) :: mesh  !< mesh partition
    integer,                intent(in) :: np(3) !< points per direction
    integer,      optional, intent(in) :: nl(3) !< number of layers

    logical, allocatable :: mask(:,:,:)
    integer, allocatable :: node(:)

    integer :: lm, lg, np1, np12, np123
    integer :: e, i, j, k, l, m, n
    logical :: partial

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
      allocate( this % master % request(size(mesh%link))    )

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
      allocate( this % ghost % request(size(mesh%link)) )

    end associate

  end subroutine Init_ElementTransferBuffer3d_Shared

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

end module Element_Transfer_Buffer_3d
