!> summary:  Generation of mesh points based on the approximate cuboids
!> author:   Joerg Stiller
!> date:     2021/03/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_GetCuboids
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generates linear mesh points based on the approximate cuboids

  module subroutine GetCuboids(mesh, x)
    class(Mesh_3D),         intent(in)  :: mesh         !< mesh partition
    real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points

    integer :: d, e, i, j, k

    !$omp single
    allocate(x(0:1, 0:1, 0:1, mesh%n_elem, 3))
    !$omp end single

    if (mesh % n_elem < 1) return

    !$omp do
    do e = 1, mesh % n_elem
      associate(x_c => mesh % element(e) % geometry % x_c)

        do d = 1, 3
          do k = 0, 1
          do j = 0, 1
          do i = 0, 1
            x(i,j,k,e,d) = x_c(0,d) + x_c(1,d) * (2*i - 1) &
                                    + x_c(2,d) * (2*j - 1) &
                                    + x_c(3,d) * (2*k - 1)
          end do
          end do
          end do
        end do

      end associate
    end do

  end subroutine GetCuboids

  !=============================================================================

end submodule MP_GetCuboids
