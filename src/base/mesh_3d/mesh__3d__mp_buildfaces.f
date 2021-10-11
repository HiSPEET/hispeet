!> summary:  Identification of mesh faces
!> author:   Joerg Stiller
!> date:     2020/12/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildFaces
  use Quick_Sort
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of mesh faces
  !>
  !> Requires
  !>   - mesh % element % vertex % id
  !>   - mesh % element % edge % id
  !>
  !> Generates
  !>   - mesh % n_face
  !>   - mesh % face
  !>   - mesh % face % element_id
  !>   - mesh % face % element_face
  !>   - mesh % element % face % id
  !>   - mesh % element % face % normal
  !>   - mesh % element % face % rotation

  module subroutine BuildFaces(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: element_face(5, 6*size(mesh%element))
    ! encoding:
    ! - element_face(1,:) :  edge 1 \
    ! - element_face(2,:) :  edge 2  } sorted according to mesh edge IDs
    ! - element_face(3,:) :  edge 3 /
    ! - element_face(4,:) :  element ID
    ! - element_face(5,:) :  corresponding element face (1 .. 6)

    integer :: nef
    integer :: i, j, k, s
    integer :: e(4), p(4)
    integer :: f1, l1
    integer :: f2, l2

    ! build ordered list of element faces ......................................

    ! extract element faces
    j = 0
    do i = 1, size(mesh%element)
      associate(element => mesh % element(i))
        do k = 1, 6
          e = element % edge(E_FACE(:,k)) % id
          call SortFaceEdges(e, p)
          element_face(:,j+k) = [ e(p(1:3)) , i, k ]
        end do
      end associate
      j = j + 6
    end do
    nef = size(element_face, 2)

    ! sort the element faces according to edge IDs
    call SortTriplets(element_face)

    ! count mesh faces .........................................................

    i = 1
    j = 1
    k = 1
    COUNT_FACES: do
      do
        j = j + 1
        if (j > nef) exit
        if (element_face(1,j) /= element_face(1,i)) exit ! edge 1 differs
        if (element_face(2,j) /= element_face(2,i)) exit ! edge 2 differs
        if (element_face(3,j) /= element_face(3,i)) exit ! edge 3 differs
      end do
      if (j > nef) exit COUNT_FACES
      i = j
      k = k + 1
    end do COUNT_FACES
    mesh % n_face = k

    ! build faces and update elements ..........................................

    allocate(mesh % face( mesh % n_face ))

    i = 1
    j = 1
    k = 1
    BUILD_FACES: do
      do
        j = j + 1
        if (j > nef) exit
        if (element_face(1,j) /= element_face(1,i)) exit ! edge 1 differs
        if (element_face(2,j) /= element_face(2,i)) exit ! edge 2 differs
        if (element_face(3,j) /= element_face(3,i)) exit ! edge 3 differs
      end do

      ! first adjacent element generates the face
      l1 = element_face(4,i) ! mesh element ID
      f1 = element_face(5,i) ! mesh element face
      select case(f1)
      case(1:3)
        s = 2   ! side 2: element in positive normal direction
      case default
        s = 1   ! side 1: element in negative normal direction
      end select
      mesh % face(k)     % element(s) % id       = l1
      mesh % face(k)     % element(s) % face     = f1
      mesh % element(l1) % face(f1)   % id       = k
      mesh % element(l1) % face(f1)   % normal   = 1
      mesh % element(l1) % face(f1)   % rotation = 0

      ! second adjacent element, if any
      if (j - i == 2) then

        ! identify second element
        l2 = element_face(4,i+1) ! mesh element ID
        f2 = element_face(5,i+1) ! mesh element face
        s  = 3 - s               ! side where the element is located

        ! update face and element data
        mesh % face(k)     % element(s) % id   = l2
        mesh % face(k)     % element(s) % face = f2
        mesh % element(l2) % face(f2)   % id   = k
        call ElementFaceOrientation(                                      &
               fv       = mesh % element(l1) % vertex(V_FACE(:,f1)) % id, &
               ev       = mesh % element(l2) % vertex(V_FACE(:,f2)) % id, &
               normal   = mesh % element(l2) % face(f2) % normal,         &
               rotation = mesh % element(l2) % face(f2) % rotation        )

      end if

      if (j > nef) exit BUILD_FACES
      i = j
      k = k + 1
    end do BUILD_FACES

  end subroutine BuildFaces

  !-----------------------------------------------------------------------------
  !> Sort face edges according to ascending ID

  pure subroutine SortFaceEdges(e, p)
    integer, intent(in)  :: e(4) !< edge IDs
    integer, intent(out) :: p(4) !< permutation: e(p(i)) ≤ e(p(j)) if i < j

    p = [1,2,3,4]
    if (e(p(1)) > e(p(2))) p([1,2]) = p([2,1])
    if (e(p(3)) > e(p(4))) p([3,4]) = p([4,3])
    if (e(p(1)) > e(p(3))) p([1,3]) = p([3,1])
    if (e(p(2)) > e(p(4))) p([2,4]) = p([4,2])
    if (e(p(2)) > e(p(3))) p([2,3]) = p([3,2])

  end subroutine SortFaceEdges

  !-----------------------------------------------------------------------------
  !> Identification of element face orientation

  pure subroutine ElementFaceOrientation(fv, ev, normal, rotation)
    integer,      intent(in)  :: fv(4)    !< mesh face vertices
    integer,      intent(in)  :: ev(4)    !< corresponding mesh element vertices
    integer(IXS), intent(out) :: normal   !< normal orientation against mesh face
    integer(IXS), intent(out) :: rotation !< 1/4-rotation to align with mesh face

    ! identify normal orientation and rotation .................................

    if (fv(1) == ev(1)) then
      if (fv(2) == ev(2)) then
        normal   =  1
        rotation =  0
      else
        normal   = -1
        rotation =  3
      end if

    else if (fv(1) == ev(2)) then
      if (fv(2) == ev(4)) then
        normal   =  1
        rotation =  1
      else
        normal   = -1
        rotation =  2
      end if

    else if (fv(1) == ev(4)) then
      if (fv(2) == ev(3)) then
        normal   =  1
        rotation =  2
      else
        normal   = -1
        rotation =  1
      end if

    else ! fv(1) == ev(3)
      if (fv(2) == ev(1)) then
        normal   =  1
        rotation =  3
      else
        normal   = -1
        rotation =  0
      end if

    end if

  end subroutine ElementFaceOrientation

  !=============================================================================

end submodule MP_BuildFaces
