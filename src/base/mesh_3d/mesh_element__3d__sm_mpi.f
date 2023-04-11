submodule(Mesh_Element__3D) SM_MPI

  !-----------------------------------------------------------------------------
  !> MPI datatype for mesh elements

  type(MPI_Datatype) :: MPI_Element = MPI_DATATYPE_NULL

  !-----------------------------------------------------------------------------
  ! MPI datatypes element components and auxiliary data

  type(MPI_Datatype) :: MPI_ElementVertex     = MPI_DATATYPE_NULL
  type(MPI_Datatype) :: MPI_ElementEdge       = MPI_DATATYPE_NULL
  type(MPI_Datatype) :: MPI_ElementFace       = MPI_DATATYPE_NULL
  type(MPI_Datatype) :: MPI_ElementAdaptation = MPI_DATATYPE_NULL
  type(MPI_Datatype) :: MPI_ElementGeometry   = MPI_DATATYPE_NULL

  logical :: initialized = .false.

contains

  !-----------------------------------------------------------------------------
  !> Get MPI data type for mesh elements
  !>
  !> The datatype allows to transfer the static components of elements, e.g.
  !>
  !>     type(MeshElement_3D) :: element(n)
  !>     type(MPI_Datatype)   :: MPI_MeshElement_3D
  !>     ...
  !>     call Get_MPI_MeshElement_3D(MPI_MeshElement_3D)
  !>     call MPI_Send(element(1)%id, n, MPI_MeshElement_3D, dest, tag, comm)

  module subroutine Get_MPI_MeshElement_3D( MPI_MeshElement_3D )
    type(MPI_Datatype), intent(out) :: MPI_MeshElement_3D

    !$omp master

    if (.not. initialized) then
      call Init_MPI_Element()
    end if

    MPI_MeshElement_3D = MPI_Element

    !$omp end master

  end subroutine Get_MPI_MeshElement_3D

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElement

  subroutine Init_MPI_Element()

    integer, parameter :: N = 9 ! number of static components

    integer(MPI_ADDRESS_KIND) :: extent  ! extent
    integer(MPI_ADDRESS_KIND) :: lb = 0  ! lower bound
    integer(MPI_ADDRESS_KIND) :: addr1   ! address of element(1)
    integer(MPI_ADDRESS_KIND) :: addr2   ! address of element(2)
    integer(MPI_ADDRESS_KIND) :: addr(N) ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N) ! displacements ...
    integer                   :: blen(N) ! block lengths ...
    type(MPI_Datatype)        :: typ(N)  ! MPI datatypes ...

    type(MPI_Datatype)   :: MPI_MeshElement_ ! provisional datatype
    type(MeshElement_3D) :: element(2)

    integer :: i

    ! preliminaries ............................................................

    ! addresses of elements(1:2)
    call MPI_Get_address(element(1) % id, addr1)
    call MPI_Get_address(element(2) % id, addr2)

    ! extent and lower bound
    extent = MPI_Aint_diff(addr2, addr1)

    ! initialize component datatypes
    call Init_MPI_ElementVertex()
    call Init_MPI_ElementEdge()
    call Init_MPI_ElementFace()
    call Init_MPI_ElementAdaptation()
    call Init_MPI_ElementGeometry()

    ! components ...............................................................

    ! adresses
    call MPI_Get_address(element(1) % id              , addr(1))
    call MPI_Get_address(element(1) % cluster_id      , addr(2))
    call MPI_Get_address(element(1) % cluster_oct     , addr(3))
    call MPI_Get_address(element(1) % frozen          , addr(4))
    call MPI_Get_address(element(1) % vertex (1) % id , addr(5))
    call MPI_Get_address(element(1) % edge   (1) % id , addr(6))
    call MPI_Get_address(element(1) % face   (1) % id , addr(7))
    call MPI_Get_address(element(1) % adaptation      , addr(8))
    call MPI_Get_address(element(1) % geometry        , addr(9))

   ! types and block lengths
    typ(1) = MPI_INTEGER           ;   blen(1) =  1
    typ(2) = MPI_INTEGER           ;   blen(2) =  1
    typ(3) = MPI_INTEGER           ;   blen(3) =  1
    typ(4) = MPI_LOGICAL           ;   blen(4) =  1
    typ(5) = MPI_ElementVertex     ;   blen(5) =  8
    typ(6) = MPI_ElementEdge       ;   blen(6) = 12
    typ(7) = MPI_ElementFace       ;   blen(7) =  6
    typ(8) = MPI_ElementAdaptation ;   blen(8) =  1
    typ(9) = MPI_ElementGeometry   ;   blen(9) =  1

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr1)
    end do

    ! resized MPI type .........................................................

    ! provisional datatype
    call MPI_Type_create_struct(N, blen, disp, typ, MPI_MeshElement_)

    ! resize and commit
    call MPI_Type_create_resized(MPI_MeshElement_, lb, extent, MPI_Element)
    call MPI_Type_commit(MPI_Element)

    ! free auxiliary datatypes
    call MPI_Type_free( MPI_MeshElement_      )
    call MPI_Type_free( MPI_ElementVertex     )
    call MPI_Type_free( MPI_ElementEdge       )
    call MPI_Type_free( MPI_ElementFace       )
    call MPI_Type_free( MPI_ElementAdaptation )
    call MPI_Type_free( MPI_ElementGeometry   )

  end subroutine Init_MPI_Element

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElementVertex

  subroutine Init_MPI_ElementVertex()

    integer, parameter :: N = 5 ! number of static components

    integer(MPI_ADDRESS_KIND) :: addr(N)     ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N)     ! displacements ...
    integer                   :: blen(N)     ! block lengths ...
    type(MPI_Datatype)        :: typ(N)      ! MPI datatypes ...

    type(MeshElementVertex_3D) :: vertex(1)

    integer :: i

    ! components ...............................................................

    ! adresses
    call MPI_Get_address(vertex(1) % id         , addr(1))
    call MPI_Get_address(vertex(1) % n_neighbor , addr(2))
    call MPI_Get_address(vertex(1) % i_neighbor , addr(3))
    call MPI_Get_address(vertex(1) % rank       , addr(4))
    call MPI_Get_address(vertex(1) % val        , addr(5))

    ! block lengths
    blen = 1

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr(1))
    end do

    ! types
    typ(1:1) = MPI_INTEGER
    typ(2:5) = MPI_INTEGER_IXS

    ! MPI datatype .............................................................

    call MPI_Type_create_struct(N, blen, disp, typ, MPI_ElementVertex)

  end subroutine Init_MPI_ElementVertex

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElementEdge

  subroutine Init_MPI_ElementEdge()

    integer, parameter :: N = 6 ! number of static components

    integer(MPI_ADDRESS_KIND) :: addr(N)     ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N)     ! displacements ...
    integer                   :: blen(N)     ! block lengths ...
    type(MPI_Datatype)        :: typ(N)      ! MPI datatypes ...

    type(MeshElementEdge_3D) :: edge(1)

    integer :: i

    ! components ...............................................................

    ! adresses
    call MPI_Get_address(edge(1) % id          , addr(1))
    call MPI_Get_address(edge(1) % orientation , addr(2))
    call MPI_Get_address(edge(1) % n_neighbor  , addr(3))
    call MPI_Get_address(edge(1) % i_neighbor  , addr(4))
    call MPI_Get_address(edge(1) % rank        , addr(5))
    call MPI_Get_address(edge(1) % val         , addr(6))

    ! block lengths
    blen = 1

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr(1))
    end do

    ! types
    typ(1:1) = MPI_INTEGER
    typ(2:6) = MPI_INTEGER_IXS

    ! MPI datatype .............................................................

    call MPI_Type_create_struct(N, blen, disp, typ, MPI_ElementEdge)

  end subroutine Init_MPI_ElementEdge

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElementFace

  subroutine Init_MPI_ElementFace()

    integer, parameter :: N = 8 ! number of static components

    integer(MPI_ADDRESS_KIND) :: addr(N)     ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N)     ! displacements ...
    integer                   :: blen(N)     ! block lengths ...
    type(MPI_Datatype)        :: typ(N)      ! MPI datatypes ...

    type(MeshElementFace_3D) :: face(1)

    integer :: i

    ! components ...............................................................

    ! adresses
    call MPI_Get_address(face(1) % id         , addr(1))
    call MPI_Get_address(face(1) % boundary   , addr(2))
    call MPI_Get_address(face(1) % normal     , addr(3))
    call MPI_Get_address(face(1) % rotation   , addr(4))
    call MPI_Get_address(face(1) % n_neighbor , addr(5))
    call MPI_Get_address(face(1) % i_neighbor , addr(6))
    call MPI_Get_address(face(1) % rank       , addr(7))
    call MPI_Get_address(face(1) % val        , addr(8))

    ! block lengths
    blen = 1

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr(1))
    end do

    ! types
    typ(1:2) = MPI_INTEGER
    typ(3:8) = MPI_INTEGER_IXS

    ! MPI datatype .............................................................

    call MPI_Type_create_struct(N, blen, disp, typ, MPI_ElementFace)

  end subroutine Init_MPI_ElementFace

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElementAdaptation

  subroutine Init_MPI_ElementAdaptation()

    integer, parameter :: N = 7 ! number of static components

    integer(MPI_ADDRESS_KIND) :: addr0       ! address of element(1)
    integer(MPI_ADDRESS_KIND) :: addr(N)     ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N)     ! displacements ...
    integer                   :: blen(N)     ! block lengths ...
    type(MPI_Datatype)        :: typ(N)      ! MPI datatypes ...

    type(MeshElementAdaptation_3D) :: adaptation

    integer :: i

    ! components ...............................................................

    ! adresse of type
    call MPI_Get_address(adaptation, addr0)

    ! adresses of components
    call MPI_Get_address(adaptation % parent_proc, addr(1))
    call MPI_Get_address(adaptation % parent_id  , addr(2))
    call MPI_Get_address(adaptation % refinement , addr(3))
    call MPI_Get_address(adaptation % sublevels  , addr(4))
    call MPI_Get_address(adaptation % child_proc , addr(5))
    call MPI_Get_address(adaptation % child_id   , addr(6))
    call MPI_Get_address(adaptation % mark       , addr(7))

    ! block lengths
    blen = 1

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr0)
    end do

    ! types
    typ = MPI_INTEGER

    ! MPI datatype .............................................................

    call MPI_Type_create_struct(N, blen, disp, typ, MPI_ElementAdaptation)

  end subroutine Init_MPI_ElementAdaptation

  !-----------------------------------------------------------------------------
  !> Initialization of MPI_MeshElementGeometry

  subroutine Init_MPI_ElementGeometry()

    integer, parameter :: N = 3 ! number of static components

    integer(MPI_ADDRESS_KIND) :: addr0       ! address of element(1)
    integer(MPI_ADDRESS_KIND) :: addr(N)     ! addresses of static components
    integer(MPI_ADDRESS_KIND) :: disp(N)     ! displacements ...
    integer                   :: blen(N)     ! block lengths ...
    type(MPI_Datatype)        :: typ(N)      ! MPI datatypes ...

    type(MeshElementGeometry_3D) :: geometry

    integer :: i

    ! components ...............................................................

    ! adresse of type
    call MPI_Get_address(geometry, addr0)

    ! adresses
    call MPI_Get_address(geometry % po   , addr(1))
    call MPI_Get_address(geometry % x_c  , addr(2))
    call MPI_Get_address(geometry % dx_m , addr(3))

    ! block lengths
    blen = [ 1, 12, 6 ]

    ! displacements
    do i = 1, N
      disp(i) = MPI_Aint_diff(addr(i), addr0)
    end do

    ! types
    typ(1:1) = MPI_INTEGER
    typ(2:3) = MPI_REAL_RNP

    ! MPI datatype .............................................................

    call MPI_Type_create_struct(N, blen, disp, typ, MPI_ElementGeometry)

  end subroutine Init_MPI_ElementGeometry

  !=============================================================================

end submodule SM_MPI
