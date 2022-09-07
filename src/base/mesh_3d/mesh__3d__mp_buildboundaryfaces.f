!> summary:  Identification of boundary faces
!> author:   Joerg Stiller
!> date:     2022/07/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildBoundaryFaces
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of boundary faces
  !>
  !> Requires
  !>   - `mesh % n_elem`
  !>   - `mesh % n_boundary`
  !>   - `mesh % boundary` without faces
  !>   - `mesh % element % face % boundary`
  !>
  !> Generates
  !>   - `mesh % boundary % n_face`
  !>   - `mesh % boundary % face`

  module subroutine BuildBoundaryFaces(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: f_bound(mesh % n_bound)
    integer :: b, e, i

    ! count local boundary faces ...............................................

    f_bound = 0

    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)
        do i = 1, 6
          b = face(i) % boundary
          if (b > 0) then
            f_bound(b) = f_bound(b) + 1
          end if
        end do
      end associate
    end do

    ! create local boundary faces ..............................................

    do b = 1, mesh % n_bound
      mesh % boundary(b) % n_face = f_bound(b)
      allocate(mesh % boundary(b) % face( f_bound(b) ))
    end do

    f_bound = 0

    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)
        do i = 1, 6
          b = face(i) % boundary
          if (b > 0) then
            f_bound(b) = f_bound(b) + 1
            mesh % boundary(b) % face(f_bound(b)) % element_id   = e
            mesh % boundary(b) % face(f_bound(b)) % element_face = i
          end if
        end do
      end associate
    end do

  end subroutine BuildBoundaryFaces

  !=============================================================================

end submodule MP_BuildBoundaryFaces
