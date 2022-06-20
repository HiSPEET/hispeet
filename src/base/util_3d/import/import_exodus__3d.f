!> summary:  Generate a generic 3d mesh from an ExodusII data base
!> author:   Susanne Stimpert, Joerg Stiller
!> date:     2016/05/29, revised 2020/12/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> This version is restricted to hexahedral elements.
!===============================================================================

module Import_Exodus__3D
  use C_Binding
  use Execution_Control
  use Exodus_Binding
  use Generic_Mesh__3D
  implicit none
  private

  public :: ImportExodus3d

contains

  !-----------------------------------------------------------------------------
  !> Reads the mesh from file and converts it to a GenericMesh_3D object

  subroutine ImportExodus3d(file, mesh)
    character(len=*),      intent(in)  :: file !< Exodus file
    class(GenericMesh_3D), intent(out) :: mesh !< output mesh

    ! Exodus interface .........................................................

    integer(C_INT) :: exo_id  ! file handle
    integer(C_INT) :: comp_ws ! application word size
    integer(C_INT) :: io_ws   ! IO word size
    real(C_FLOAT)  :: version ! EXODUS II database version number
    integer(C_INT) :: err     ! error code

    ! dimensions
    integer(C_INT) :: num_dim, num_nodes, num_elem, num_elem_blk
    integer(C_INT) :: num_node_sets, num_side_sets

    character(C_CHAR) :: title(EX_MAX_LINE_LENGTH+1) !

    real(C_DOUBLE), allocatable :: xn(:), yn(:), zn(:)
    integer(C_INT), allocatable :: elem_blk_ids(:)

    character(C_CHAR) :: elem_type_this_blk(EX_MAX_STR_LENGTH+1)
    integer(C_INT)    :: num_elem_this_blk, num_nodes_per_elem, num_attr

    character(len=EX_MAX_STR_LENGTH) :: elem_type
    logical :: linear    = .false.  ! switch for tri-linear hexahedra
    logical :: quadratic = .false.  ! switch for tri-quadratic hexahedra

    integer(C_INT), allocatable :: connect(:,:) ! element connectivity lists

    integer(C_INT) :: num_side_in_set, num_dist_fact_in_set
    integer(C_INT), allocatable :: side_set_ids(:)
    integer(C_INT), allocatable :: side_set_elem_list(:)
    integer(C_INT), allocatable :: side_set_side_list(:)

    character(C_CHAR) :: side_set_name(EX_MAX_STR_LENGTH+1)

    ! Auxiliary variables ........................................................

    integer, allocatable :: vert(:) ! maps node to vertex IDs
    integer :: num_vert

    integer :: i, k, l
    logical :: failed

    ! Open Exodus file and get parameters ........................................

    ! set application word size to size of double
    comp_ws = c_sizeof(1.0_C_DOUBLE)

    ! set io_ws to 0, i.e., adopt the IO word size as used in the Exodus file
    io_ws = 0

    ! open Exodus file
    write(*,'(/,A,/)') 'Opening Exodus II file '//trim(file)
    exo_id = ExOpen(CString(file), EX_READ, comp_ws, io_ws, version)

    ! initialization parameters
    call ExGetInit( exo_id, title, num_dim, num_nodes, num_elem,    &
                    num_elem_blk, num_node_sets, num_side_sets, err )

    write(*,'(2X,2A,/)') 'title: ', trim(FString(title))
    write(*,'(2X,A,I0)') 'number of nodes      = ', num_nodes
    write(*,'(2X,A,I0)') 'number of elements   = ', num_elem
    write(*,'(2X,A,I0)') 'number of blocks     = ', num_elem_blk
    write(*,'(2X,A,I0)') 'number of boundaries = ', num_side_sets
    write(*,*)

    if (num_elem < 1) then
      allocate(mesh%vertex(0))
      allocate(mesh%element(0))
      allocate(mesh%boundary(0))
      call Warning('ImportExodusMesh','Returning empty mesh')
      return
    end if

    ! element block IDs
    allocate(elem_blk_ids(num_elem_blk))
    call ExGetElemBlkIds(exo_id, elem_blk_ids, err)

    ! read properties of first element block
    call ExGetElemBlock( exo_id, elem_blk_ids(1), elem_type_this_blk,         &
                         num_elem_this_blk, num_nodes_per_elem, num_attr, err )

    ! identify the element type
    elem_type = FString(elem_type_this_blk)
    if (elem_type(1:3) == 'HEX') then
      linear    = num_nodes_per_elem == 8
      quadratic = num_nodes_per_elem == 27
    end if
    if (linear) then
      write(*,'(2X,A,/)') 'importing tri-linear hexahedral elements ...'
    else if (quadratic) then
      write(*,'(2X,A,/)') 'importing tri-quadratic hexahedral elements ...'
    else
      call Error('ImportExodus', 'Only HEX8 and HEX27 elements are supported')
    end if

    ! read node coordinates ....................................................

    allocate(xn(num_nodes), yn(num_nodes), zn(num_nodes))
    call ExGetCoord(exo_id, xn, yn, zn, err)

    ! read element connectivity ................................................

    if (linear) then
      allocate(connect( 8, num_elem), source = -1_C_INT)
    else ! quadratic
      allocate(connect(27, num_elem), source = -1_C_INT)
    end if

    l = 1
    do i = 1, num_elem_blk

      call ExGetElemBlock( exo_id, elem_blk_ids(i), elem_type_this_blk,     &
                           num_elem_this_blk, num_nodes_per_elem, num_attr, &
                           err                                              )

      elem_type = FString(elem_type_this_blk)
      if (linear) then
        failed = elem_type(1:3) /= 'HEX' .or. num_nodes_per_elem /= 8
      else ! quadratic
        failed = elem_type(1:3) /= 'HEX' .or. num_nodes_per_elem /= 27
      end if

      call ExGetElemConn(exo_id, elem_blk_ids(i), connect(:,l:), err)

      l = l + num_elem_this_blk

    end do

    ! identify vertices .......................................................

    allocate(vert(num_nodes), source = 0)

    ! mark element vertex nodes
    do i = 1, num_elem
      vert( connect(1:8,i) ) = 1
    end do

    ! identify mesh vertices
    k = 0
    do i = 1, num_nodes
      if (vert(i) == 0) cycle
      k = k + 1
      vert(i) = k
    end do
    num_vert = k

    write(*,'(2X,3(A,1X,I0,1X))') &
      'extracted', num_vert, 'mesh vertices from', num_nodes, 'nodes'
    write(*,'(/,2X,A)') 'bounding box:'
    write(*,'(4X,A,2(ES10.3,2X,A))') 'x = ', minval(xn), ' ... ', maxval(xn)
    write(*,'(4X,A,2(ES10.3,2X,A))') 'y = ', minval(yn), ' ... ', maxval(yn)
    write(*,'(4X,A,2(ES10.3,2X,A))') 'z = ', minval(zn), ' ... ', maxval(zn)
    write(*,*)

    ! create generic mesh vertices .............................................

    allocate (mesh%vertex(num_vert))
    do i = 1, num_nodes
      k = vert(i)
      if (k == 0) cycle
      mesh%vertex(k)%id = k
      mesh%vertex(k)%x  = [ xn(i), yn(i), zn(i) ]
    end do

    ! create mesh elements .....................................................

    ! set rotational numbering
    mesh%numbering = ROTATIONAL_NUMBERING

    allocate(mesh%element(num_elem))

    do l = 1, num_elem

      ! ID and type
      mesh%element(l)%id  = l
      mesh%element(l)%typ = HEXAHEDRAL_ELEMENT

      ! vertices
      do k = 1, 8
        mesh%element(l)%vertex(k) = vert( connect(k,l) )
      end do

      ! mapping to physical element domain
      if (linear) then
        mesh%element(l)%order = 1
        mesh%element(l)%basis = GAUSS_LOBATTO_BASIS
        allocate(mesh%element(l)%x(8,3))
        associate(x => mesh%element(l)%x, node => connect(:,l))
          x(1, 1:3) = [ xn( node(1) ), yn( node(1) ), zn( node(1) ) ]
          x(2, 1:3) = [ xn( node(2) ), yn( node(2) ), zn( node(2) ) ]
          x(3, 1:3) = [ xn( node(4) ), yn( node(4) ), zn( node(4) ) ]
          x(4, 1:3) = [ xn( node(3) ), yn( node(3) ), zn( node(3) ) ]
          x(5, 1:3) = [ xn( node(5) ), yn( node(5) ), zn( node(5) ) ]
          x(6, 1:3) = [ xn( node(6) ), yn( node(6) ), zn( node(6) ) ]
          x(7, 1:3) = [ xn( node(8) ), yn( node(8) ), zn( node(8) ) ]
          x(8, 1:3) = [ xn( node(7) ), yn( node(7) ), zn( node(7) ) ]
        end associate

      else if (quadratic) then

        mesh%element(l)%order = 2
        mesh%element(l)%basis = GAUSS_LOBATTO_BASIS
        allocate(mesh%element(l)%x(27,3))
        associate(x => mesh%element(l)%x, node => connect(:,l))
          x( 1, 1:3) = [ xn( node( 1) ), yn( node( 1) ), zn( node( 1) ) ]
          x( 2, 1:3) = [ xn( node( 9) ), yn( node( 9) ), zn( node( 9) ) ]
          x( 3, 1:3) = [ xn( node( 2) ), yn( node( 2) ), zn( node( 2) ) ]
          x( 4, 1:3) = [ xn( node(12) ), yn( node(12) ), zn( node(12) ) ]
          x( 5, 1:3) = [ xn( node(22) ), yn( node(22) ), zn( node(22) ) ]
          x( 6, 1:3) = [ xn( node(10) ), yn( node(10) ), zn( node(10) ) ]
          x( 7, 1:3) = [ xn( node( 4) ), yn( node( 4) ), zn( node( 4) ) ]
          x( 8, 1:3) = [ xn( node(11) ), yn( node(11) ), zn( node(11) ) ]
          x( 9, 1:3) = [ xn( node( 3) ), yn( node( 3) ), zn( node( 3) ) ]
          x(10, 1:3) = [ xn( node(13) ), yn( node(13) ), zn( node(13) ) ]
          x(11, 1:3) = [ xn( node(26) ), yn( node(26) ), zn( node(26) ) ]
          x(12, 1:3) = [ xn( node(14) ), yn( node(14) ), zn( node(14) ) ]
          x(13, 1:3) = [ xn( node(24) ), yn( node(24) ), zn( node(24) ) ]
          x(14, 1:3) = [ xn( node(21) ), yn( node(21) ), zn( node(21) ) ]
          x(15, 1:3) = [ xn( node(25) ), yn( node(25) ), zn( node(25) ) ]
          x(16, 1:3) = [ xn( node(16) ), yn( node(16) ), zn( node(16) ) ]
          x(17, 1:3) = [ xn( node(27) ), yn( node(27) ), zn( node(27) ) ]
          x(18, 1:3) = [ xn( node(15) ), yn( node(15) ), zn( node(15) ) ]
          x(19, 1:3) = [ xn( node( 5) ), yn( node( 5) ), zn( node( 5) ) ]
          x(20, 1:3) = [ xn( node(17) ), yn( node(17) ), zn( node(17) ) ]
          x(21, 1:3) = [ xn( node( 6) ), yn( node( 6) ), zn( node( 6) ) ]
          x(22, 1:3) = [ xn( node(20) ), yn( node(20) ), zn( node(20) ) ]
          x(23, 1:3) = [ xn( node(23) ), yn( node(23) ), zn( node(23) ) ]
          x(24, 1:3) = [ xn( node(18) ), yn( node(18) ), zn( node(18) ) ]
          x(25, 1:3) = [ xn( node( 8) ), yn( node( 8) ), zn( node( 8) ) ]
          x(26, 1:3) = [ xn( node(19) ), yn( node(19) ), zn( node(19) ) ]
          x(27, 1:3) = [ xn( node( 7) ), yn( node( 7) ), zn( node( 7) ) ]
        end associate

      end if

    end do

    ! read side sets and create boundaries .....................................

    allocate(mesh%boundary(num_side_sets))

    ! side set IDs
    allocate(side_set_ids(num_side_sets))
    call ExGetSideSetIds(exo_id, side_set_ids, err)

    do l = 1, num_side_sets

      ! identifier
      mesh%boundary(l)%id = l

      ! name
      call ExGetName(exo_id, EX_SIDE_SET, side_set_ids(l), side_set_name, err)
      mesh%boundary(l)%name = FString(side_set_name)

      ! number of faces = num_side_in_set
      call ExGetSideSetParam( exo_id, side_set_ids(l), num_side_in_set, &
                              num_dist_fact_in_set, err                 )

      ! print status info
      write(*,'(2X,A,I0,1X,A,1X,A,I0,1X,A)') &
        'importing boundary ', l, &
        'named "' // trim(mesh%boundary(l)%name) // '"', &
        'with ', num_side_in_set, 'faces'

      ! element and element face lists
      allocate(side_set_elem_list( num_side_in_set ))
      allocate(side_set_side_list( num_side_in_set ))

      ! get element and face ids
      call ExGetSideSet( exo_id, side_set_ids(l), side_set_elem_list, &
                         side_set_side_list, err                      )

      ! create boundary faces
      allocate(mesh%boundary(l)%face(num_side_in_set))
      do i = 1, num_side_in_set
        mesh%boundary(l)%face(i)%element_id   = side_set_elem_list(i)
        mesh%boundary(l)%face(i)%element_face = side_set_side_list(i)
      end do

      ! release workspace
      deallocate(side_set_elem_list, side_set_side_list)

    end do

    write(*,*)

    ! cleanup ..................................................................

    call ExClose(exo_id, err)

  end subroutine ImportExodus3d

  !=============================================================================

end module Import_Exodus__3D
