!> summary:  3D ML INS time averaging
!> author:   Joerg Stiller
!> date:     2025/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__INS__Time_Averaging__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  use Execution_Control
  use ML__Mesh_Variable__3D

  implicit none
  private

  public :: ML_INS_TimeAveraging_3D

contains

  !-----------------------------------------------------------------------------
  !> Multilevel INS time averaging
  !>
  !> Computes
  !>   - temporal averages of the flow variables contained in `u`
  !>   - temporal averages of the squared flow variables
  !>   - the upper triangular part of the Reynolds stress tensor

  subroutine ML_INS_TimeAveraging_3D(u, a, n)
    class(ML_MeshVariable_3D), intent(in)    :: u  !< flow variables
    class(ML_MeshVariable_3D), intent(inout) :: a  !< averaged quantities
    integer,                   intent(inout) :: n  !< number of samples

    real(RNP) :: w_a, w_s
    integer   :: c, e, l

    n = n + 1
    if (n > 1) then
      w_a = real(n-1, RNP) / n ! weight of current average
      w_s = ONE - w_a          ! weight of sample
    else
      w_a = 0
      w_s = 1
    end if

    do l = 1, size(u % level)
      associate( mesh => u % level(l) % mesh &
               , a_l  => a % level(l) % val  &
               , u_l  => u % level(l) % val  &
               , nc   => u % level(l) % nc   )

        !$omp do
        do e = 1, mesh % n_elem

          do c = 1, nc
            a_l(:,:,:,e,c   ) = w_a * a_l(:,:,:,e,c   ) + w_s * u_l(:,:,:,e,c)
            a_l(:,:,:,e,c+nc) = w_a * a_l(:,:,:,e,c+nc) + w_s * u_l(:,:,:,e,c)**2
          end do

          a_l(:,:,:,e,2*nc+1) = w_a * a_l(:,:,:,e,2*nc+1) &
                              + w_s * u_l(:,:,:,e,1) * u_l(:,:,:,e,2)
          a_l(:,:,:,e,2*nc+2) = w_a * a_l(:,:,:,e,2*nc+2) &
                              + w_s * u_l(:,:,:,e,1) * u_l(:,:,:,e,3)
          a_l(:,:,:,e,2*nc+3) = w_a * a_l(:,:,:,e,2*nc+3) &
                              + w_s * u_l(:,:,:,e,2) * u_l(:,:,:,e,3)

        end do

      end associate
    end do

  end subroutine ML_INS_TimeAveraging_3D

  !=============================================================================

end module ML__INS__Time_Averaging__3D
