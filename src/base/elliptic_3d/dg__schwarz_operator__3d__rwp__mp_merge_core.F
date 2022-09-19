  !-----------------------------------------------------------------------------
  !>  Merge of subdomain core contributions into mesh variable

  module subroutine Merge_Core(this, mesh, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
    real(RWP),      intent(inout) :: us(:,:,:,:)   !< mesh variable

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
          u(i,j,k,e) = u(i,j,k,e) + us(c0 + i, c0 + j, c0 + k, e)
        end do
        end do
        end do
      end do

    end associate

  end subroutine Merge_Core
