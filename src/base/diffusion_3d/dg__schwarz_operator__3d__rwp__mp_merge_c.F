  !-----------------------------------------------------------------------------
  !>  Merge of subdomain core contributions into mesh variable

  module subroutine Merge_C(this, mesh, vs, v)
    class(DG_SchwarzOperator_3D_RWP), intent(in) :: this
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< subdomain solutions
    real(RWP),      intent(inout) :: vs(:,:,:,:)   !< mesh variable

    integer :: e, i, j, k
    integer :: c0

    associate(po => this%po, no => this%no)

      ! index of first core point within subdomain
      c0 = no + 1

      !$omp do
      do e = 1, mesh % n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          v(i,j,k,e) = v(i,j,k,e) + vs(c0 + i, c0 + j, c0 + k, e)
        end do
        end do
        end do
      end do

    end associate

  end subroutine Merge_C
