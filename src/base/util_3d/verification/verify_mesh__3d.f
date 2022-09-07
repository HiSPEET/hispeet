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

    logical :: con_passed

    call ElementConnectivityTest(mesh, con_passed)

    passed = con_passed

  end subroutine VerifyMesh_3D

  !-----------------------------------------------------------------------------
  !> Test of local element connectivity

  subroutine ElementConnectivityTest(mesh, passed)
    type(Mesh_3D), intent(in)  :: mesh   !< mesh partition
    logical,       intent(out) :: passed !< test result

    integer :: c, e, i, j, k, l, n

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
              passed = FaceMatch(mesh%element(l), c, e)
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
              passed = EdgeMatch(mesh%element(l), c, e)
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
              passed = VertexMatch(mesh%element(l), c, e)
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

  logical function FaceMatch(element, ef, ml) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ef !< coupled element face
    integer, intent(in) :: ml !< mesh element ID to match

    integer :: i, j, n

    match = .false.

    n = element % face(ef) % n_neighbor
    i = element % face(ef) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (match) exit
    end do

  end function FaceMatch

  !-----------------------------------------------------------------------------
  !> Edge matching test

  logical function EdgeMatch(element, ee, ml) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ee !< coupled element edge
    integer, intent(in) :: ml !< mesh element ID to match

    integer :: i, j, n

    match = .false.

    n = element % edge(ee) % n_neighbor
    i = element % edge(ee) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (match) exit
    end do

  end function EdgeMatch

  !-----------------------------------------------------------------------------
  !> Vertex matching test

  logical function VertexMatch(element, ev, ml) result(match)
    class(MeshElement_3D), intent(in) :: element !< given element
    integer, intent(in) :: ev !< coupled element vertex
    integer, intent(in) :: ml !< mesh element ID to match

    integer :: i, j, n

    match = .false.

    n = element % vertex(ev) % n_neighbor
    i = element % vertex(ev) % i_neighbor
    do j = i, i+n-1
      match = element % neighbor(j) % id == ml
      if (match) exit
    end do

  end function VertexMatch

  !-----------------------------------------------------------------------------
  !> Averaging check

  !=============================================================================

end module Verify_Mesh__3D

