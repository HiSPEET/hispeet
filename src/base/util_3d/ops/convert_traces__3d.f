!> summary:  Conversion of interior traces to exterior traces
!> author:   Joerg Stiller
!> date:     2022/01/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Convert_Traces__3D
  use Kind_Parameters, only: RNP
  use Mesh__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Converts interior traces u⁻ to exterior traces u⁺
  !>
  !> This routine converts the interior traces `u⁻ = um` into exterior traces
  !> `u⁺ = up`. Both are aligned with the element faces. In order to transfer
  !> the traces of remote elements, they must be present in the corresponding
  !> ghost entries of `um`.

  subroutine ConvertTraces(mesh, um, up)
    class(Mesh_3D), intent(in) :: mesh
    real(RNP), contiguous, intent(in)  :: um(:,:,:,:,:) !< u⁻
    real(RNP), contiguous, intent(out) :: up(:,:,:,:,:) !< u⁺

    integer :: c, e, i, f, l, m, nc

    nc = size(um,5)

    !$omp do
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))
        do f = 1, 6
          i = element % face(f) % i_neighbor
          if (i > 0) then ! adopt trace from neighbor
            l = element % neighbor(i) % id
            m = element % neighbor(i) % component
            do c = 1, nc
              call element % AlignFromNeighborFace &
                                 (f, i, um(:,:,m,l,c), up(:,:,f,e,c))
            end do
          else ! no neighbor
            do c = 1, nc
              up(:,:,f,e,c) = um(:,:,f,e,c)
            end do
          end if
        end do
      end associate
    end do

  end subroutine ConvertTraces

  !=============================================================================

end module Convert_Traces__3D
