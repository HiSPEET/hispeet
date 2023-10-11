!> summary:  Verification of 3D mesh data
!> author:   Joerg Stiller
!> date:     2021/02/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Verify_Mesh__3D

  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D

  implicit none
  private

  public :: VerifyMesh_3D

contains

  !-----------------------------------------------------------------------------
  !> 3D mesh data verification procedure

  subroutine VerifyMesh_3D(mesh, passed)
    type(Mesh_3D), intent(in)  :: mesh !< mesh partition
    logical,       intent(out) :: passed !< test result

    logical :: passed_connection
    logical :: passed_enclosure

    call ElementConnectivityTest (mesh, passed_connection)
    call ElementEnclosureTest    (mesh, passed_enclosure)

    passed = passed_connection .and. passed_enclosure

  end subroutine VerifyMesh_3D

  !-----------------------------------------------------------------------------
  !> Test of local element connectivity

  subroutine ElementConnectivityTest(mesh, passed)
    type(Mesh_3D), intent(in)  :: mesh   !< mesh partition
    logical,       intent(out) :: passed !< test result

    integer :: c, e, i, j, k, l, n, o

    passed = .true.

    ELEMENTS: do e = 1, mesh % n_elem
      associate( face     => mesh % element(e) % face     &
               , edge     => mesh % element(e) % edge     &
               , vertex   => mesh % element(e) % vertex   &
               , neighbor => mesh % element(e) % neighbor )

        ! faces
        do k = 1, 6
          n = face(k) % n_neighbor
          i = face(k) % i_neighbor
          do j = i, i+n-1
            l = neighbor(j) % id
            if (l > 0 .and. l <= mesh % n_elem) then
              c = ElementFaceID(neighbor(j) % component)
              o = ElementFaceID(neighbor(j) % orientation)
              passed = FaceMatch(mesh%element(l), c, e, o)
              if (.not. passed) then
                write(*,'(99G0)') '*** part ', mesh%part,': no match between ', &
                                  'element ', e, ', face ', k, ' and ',         &
                                  'element ', l, ', face ', c
              end if
            else if (l < 0 .or. l > mesh%n_elem + mesh%n_ghost) then
              passed = .false.
              write(*,'(99G0)') '*** part ', mesh%part,': invalid neighbor at ', &
                                'element ', e, ', face ', k
            end if
            if (.not. passed) exit ELEMENTS
          end do
        end do

        ! edges
        do k = 1, 12
          n = edge(k) % n_neighbor
          i = edge(k) % i_neighbor
          do j = i, i+n-1
            l = neighbor(j) % id
            if (l > 0 .and. l <= mesh % n_elem) then
              c = ElementEdgeID(neighbor(j) % component)
              o = ElementFaceID(neighbor(j) % orientation)
              passed = EdgeMatch(mesh%element(l), c, e, o)
              if (.not. passed) then
                write(*,'(99G0)') '*** part ', mesh%part,': no match between ', &
                                  'element ', e, ', edge ', k, ' and ',         &
                                  'element ', l, ', edge ', c
              end if
            else if (l < 0 .or. l > mesh%n_elem + mesh%n_ghost) then
              passed = .false.
              write(*,'(99G0)') '*** part ', mesh%part,': invalid neighbor at ', &
                                'element ', e, ', edge ', k
            end if
            if (.not. passed) exit ELEMENTS
          end do
        end do

        ! vertices
        do k = 1, 8
          n = vertex(k) % n_neighbor
          i = vertex(k) % i_neighbor
          do j = i, i+n-1
            l = neighbor(j) % id
            if (l > 0 .and. l <= mesh % n_elem) then
              c = ElementVertexID(neighbor(j) % component)
              o = ElementFaceID(neighbor(j) % orientation)
              passed = VertexMatch(mesh%element(l), c, e, o)
              if (.not. passed) then
                write(*,'(99G0)') '*** part ', mesh%part,': no match between ',  &
                                  'element ', e, ', vertex ', k, ' and ',        &
                                  'element ', l, ', vertex ', c
              end if
            else if (l < 0 .or. l > mesh%n_elem + mesh%n_ghost) then
              passed = .false.
              write(*,'(99G0)') '*** part ', mesh%part,': invalid neighbor at ', &
                                'element ', e, ', vertex ', k
            end if
            if (.not. passed) exit ELEMENTS
          end do
        end do

      end associate
    end do ELEMENTS

  end subroutine ElementConnectivityTest

  !-----------------------------------------------------------------------------
  !> Face matching test

  logical function FaceMatch(element, ef, ml, no) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ef !< coupled element face
    integer, intent(in) :: ml !< mesh element ID to match
    integer, intent(in) :: no !< neighbor orientation

    integer :: i, j, n

    match = .false.

    n = element % face(ef) % n_neighbor
    i = element % face(ef) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (.not. match) cycle
      match = element % neighbor(j) % orientation == ReverseOrientation(no)
      if (match) exit
    end do

  end function FaceMatch

  !-----------------------------------------------------------------------------
  !> Edge matching test

  logical function EdgeMatch(element, ee, ml, no) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ee !< coupled element edge
    integer, intent(in) :: ml !< mesh element ID to match
    integer, intent(in) :: no !< neighbor orientation

    integer :: i, j, n

    match = .false.

    n = element % edge(ee) % n_neighbor
    i = element % edge(ee) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (.not. match) cycle
      match = element % neighbor(j) % orientation == ReverseOrientation(no)
      if (match) exit
    end do

  end function EdgeMatch

  !-----------------------------------------------------------------------------
  !> Vertex matching test

  logical function VertexMatch(element, ev, ml, no) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ev !< coupled element vertex
    integer, intent(in) :: ml !< mesh element ID to match
    integer, intent(in) :: no !< neighbor orientation

    integer :: i, j, n

    match = .false.

    n = element % vertex(ev) % n_neighbor
    i = element % vertex(ev) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (.not. match) cycle
      match = element % neighbor(j) % orientation == ReverseOrientation(no)
      if (match) exit
    end do

  end function VertexMatch

  !-----------------------------------------------------------------------------
  ! Returns the reverse orientation

  integer function ReverseOrientation(o) result(r)
    integer, intent(in) :: o

    select case(o)

    case(12);  r = 12
    case(13);  r = 16
    case(15);  r = 15
    case(16);  r = 13

    case(21);  r = 21
    case(23);  r = 31
    case(24);  r = 51
    case(26);  r = 61

    case(31);  r = 23
    case(32);  r = 62
    case(34);  r = 56
    case(35);  r = 35

    case(42);  r = 42
    case(43);  r = 43
    case(45);  r = 45
    case(46);  r = 46

    case(51);  r = 24
    case(53);  r = 64
    case(54);  r = 54
    case(56);  r = 34

    case(61);  r = 26
    case(62);  r = 32
    case(64);  r = 53
    case(65);  r = 65

    case default
      r = -1

    end select

  end function ReverseOrientation

  !-----------------------------------------------------------------------------
  !> Enclosure check
  !>
  !> For each active element, check if there is a neighbor on all interior
  !> faces. For frozen elements, check if they have at least one neighbor.

  subroutine ElementEnclosureTest(mesh, passed)
    type(Mesh_3D), intent(in)  :: mesh   !< mesh partition
    logical,       intent(out) :: passed !< test result

    integer :: e, k

    passed = .true.

    ELEMENTS: do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        if (element % frozen) then

          if (size(element % neighbor) < 1) then
            passed = .false.
            write(*,'(99G0)') '*** part ', mesh%part, &
                              ': isolated frozen element ', e
          end if

        else

          do k = 1, 6
            if (element % face(k) % boundary   > 0) cycle
            if (element % face(k) % n_neighbor < 1) then
              passed = .false.
              write(*,'(99G0)') '*** part ', mesh%part, &
                                ': missing neighbor at element ',e,', face ',k
!### CHECK
!! write(*,'(99(G0,1X))') 'mesh%element(',e,')%adaptation%parent_id =', &
!!                         mesh%element(  e  )%adaptation%parent_id
!### CHECK
            end if
          end do

        end if

      end associate
    end do ELEMENTS

  end subroutine ElementEnclosureTest

  !=============================================================================

end module Verify_Mesh__3D

