!> summary:  Identification of element component ranks
!> author:   Joerg Stiller
!> date:     2021/02/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_IdentifyComponentRanks
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of element component ranks
  !>
  !> Requires
  !>   - mesh % element % face   % {n_neighbor, i_neighbor}
  !>   - mesh % element % edge   % {n_neighbor, i_neighbor}
  !>   - mesh % element % vertex % {n_neighbor, i_neighbor}
  !>   - mesh % element % neighbor
  !>
  !> Generates
  !>   - mesh % element % face   % {rank, val}
  !>   - mesh % element % edge   % {rank, val}
  !>   - mesh % element % vertex % {rank, val}

  module subroutine IdentifyComponentRanks(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition

    integer :: e, k

    !$omp do schedule(static)
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        do k = 1, 6
          call element % DetermineFaceRank(k, rank = element%face(k)%rank, &
                                              val  = element%face(k)%val   )
        end do

        do k = 1, 12
          call element % DetermineEdgeRank(k, rank = element%edge(k)%rank, &
                                              val  = element%edge(k)%val   )
        end do

        do k = 1, 8
          call element % DetermineVertexRank(k, rank = element%vertex(k)%rank, &
                                                val  = element%vertex(k)%val   )
        end do

      end associate
    end do

  end subroutine IdentifyComponentRanks

  !=============================================================================

end submodule MP_IdentifyComponentRanks
