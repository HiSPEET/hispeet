!> summary:  Identification of primary element components
!> author:   Joerg Stiller
!> date:     2021/02/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_IdentifyPrimaries
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of primary element components
  !>
  !> Requires
  !>   - mesh % element % face   % {n_neighbor, i_neighbor}
  !>   - mesh % element % edge   % {n_neighbor, i_neighbor}
  !>   - mesh % element % vertex % {n_neighbor, i_neighbor}
  !>   - mesh % element % neighbor
  !>
  !> Generates
  !>   - mesh % element % face   % primary
  !>   - mesh % element % edge   % primary
  !>   - mesh % element % vertex % primary

  module subroutine IdentifyPrimaries(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition

    integer :: e, i, j, k, n

    !$omp do schedule(static)
    do e = 1, mesh % n_elem
      associate( face     => mesh % element(e) % face     &
               , edge     => mesh % element(e) % edge     &
               , vertex   => mesh % element(e) % vertex   &
               , neighbor => mesh % element(e) % neighbor )

        Faces: do k = 1, 6
          n = face(k) % n_neighbor
          i = face(k) % i_neighbor
          face(k) % primary = 1_IXS
          do j = i, i+n-1
            if (e > neighbor(j) % id) then
              face(k) % primary = 0_IXS
              exit
            end if
          end do
        end do Faces

        Edges: do k = 1, 12
          n = edge(k) % n_neighbor
          i = edge(k) % i_neighbor
          edge(k) % primary = 1_IXS
          do j = i, i+n-1
            if (e > neighbor(j) % id) then
              edge(k) % primary = 0_IXS
              exit
            end if
          end do
        end do Edges

        Vertices: do k = 1, 12
          n = vertex(k) % n_neighbor
          i = vertex(k) % i_neighbor
          vertex(k) % primary = 1_IXS
          do j = i, i+n-1
            if (e > neighbor(j) % id) then
              vertex(k) % primary = 0_IXS
              exit
            end if
          end do
        end do Vertices

      end associate
    end do

  end subroutine IdentifyPrimaries

  !=============================================================================

end submodule MP_IdentifyPrimaries
