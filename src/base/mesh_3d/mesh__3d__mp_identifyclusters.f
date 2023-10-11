!> summary:  Identification of sibling element clusters
!> author:   Joerg Stiller
!> date:     2023/04/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_IdentifyClusters
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of sibling element clusters
  !>
  !> Requires
  !>   - mesh % is_root
  !>   - mesh % element % adaptation % parent_proc
  !>   - mesh % element % adaptation % parent_id
  !>
  !> Generates
  !>   - mesh % n_cluster
  !>   - mesh % element % cluster_id

  module subroutine IdentifyClusters(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    integer :: e, parent_id, parent_proc

    !$omp master
    associate(c => mesh % n_cluster)
      c = 0
      parent_id   = -1
      parent_proc = -1
      do e = 1, mesh % n_elem
        if ( mesh % element(e) % adaptation % parent_proc /= parent_proc .or. &
             mesh % element(e) % adaptation % parent_id   /= parent_id      ) &
        then
          parent_proc = mesh % element(e) % adaptation % parent_proc
          parent_id   = mesh % element(e) % adaptation % parent_id
          c = c + 1
        end if
        mesh % element(e) % cluster_id = c
      end do
    end associate
    !$omp end master
    !$omp barrier

  end subroutine IdentifyClusters

  !=============================================================================

end submodule MP_IdentifyClusters
