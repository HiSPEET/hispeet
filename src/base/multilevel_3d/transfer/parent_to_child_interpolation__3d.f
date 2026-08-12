!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Interpolation of mesh data from parent to child level
!> author:   Joerg Stiller
!> date:     2024/08/05
!===============================================================================

module Parent_To_Child_Interpolation__3D
  use Kind_Parameters
  use Array_Assignments
  use HP__Refinement_Operator__1D
  use Mesh__3D
  use Data_Exchange__3D
  use TPO__AAA__3D
  use TPO__1To8__3D
  implicit none
  private

  public :: ParentToChildInterpolation_3D

  interface ParentToChildInterpolation_3D
    module procedure ParentToChildInterpolation_A
    module procedure ParentToChildInterpolation_S
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D parent to child interpolation of array data

  subroutine ParentToChildInterpolation_A( parent, child, iop, v_p, v_c &
                                         , skip_frozen                  )
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< interpolation operator
    real(RNP), contiguous, intent(in) :: v_p(0:,0:,0:,:,:)
      !< parent data
    real(RNP), contiguous, intent(inout) :: v_c(0:,0:,0:,:,:)
      !< child data
    logical, optional, intent(in) :: skip_frozen
      !< exclude frozen elements [F]

    call ParentToChildInterpolation_X( parent, child, iop, size(v_p,5) &
                                     , v_p, v_c, skip_frozen           )

  end subroutine ParentToChildInterpolation_A

  !-----------------------------------------------------------------------------
  !> 3D parent to child interpolation of scalar data

  subroutine ParentToChildInterpolation_S( parent, child, iop, v_p, v_c &
                                         , skip_frozen                  )
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< interpolation operator
    real(RNP), contiguous, intent(in) :: v_p(0:,0:,0:,:)
      !< parent data
    real(RNP), contiguous, intent(inout) :: v_c(0:,0:,0:,:)
      !< child data
    logical, optional, intent(in) :: skip_frozen
      !< exclude frozen elements [F]

    call ParentToChildInterpolation_X( parent, child, iop, 1 &
                                     , v_p, v_c, skip_frozen )

  end subroutine ParentToChildInterpolation_S

  !-----------------------------------------------------------------------------
  !> 3D parent to child interpolation using arrays of explicit shape

  subroutine ParentToChildInterpolation_X( parent, child, iop, n_comp &
                                         , v_p, v_c, skip_frozen      )

    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< interpolation operator
    integer, intent(in) :: n_comp
      !< number of array components
    real(RNP), intent(in) :: &
      v_p(0:iop%po_c, 0:iop%po_c, 0:iop%po_c, parent% n_elem, n_comp)
      !< parent data
    real(RNP), target, intent(inout) :: &
      v_c(0:iop%po_f, 0:iop%po_f, 0:iop%po_f, child % n_elem, n_comp)
      !< child data
    logical, optional, intent(in) :: skip_frozen
      !< exclude frozen elements [F]

    type(DataExchangeMap_3D), allocatable, save :: send_map(:)
    type(DataExchangeMap_3D), allocatable, save :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable, save :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable, save :: recv_buf(:)

    real(RNP), allocatable,         save :: v_r(:,:,:,:,:)
    real(RNP), contiguous, pointer, save :: v_i(:,:,:,:,:)

    integer, save :: n_recv, n_send

    logical :: intermediate_data
    integer :: e, i, k

    ! initialization ...........................................................

    ! check intermediate data is required to handle frozen elements
    intermediate_data = parent%refinement == 's' .and. child%n_elem_frozen > 0
    if (present(skip_frozen)) then
      intermediate_data = intermediate_data .and. .not. skip_frozen
    end if

    ! parent send map and send buffer
    n_send = parent % n_child
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(parent % map_child(i), skip_frozen)
    end do

    ! child receive map and receive buffer
    n_recv = child % n_parent
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(child % map_parent(i), skip_frozen)
    end do

    ! prepare parent data received by child
    allocate(v_r(0:iop%po_c,0:iop%po_c,0:iop%po_c,child%n_cluster,n_comp))

    ! prepare intermediate data
    if (intermediate_data) then
      allocate(v_i(0:iop%po_c,0:iop%po_c,0:iop%po_c,8*child%n_cluster,n_comp))
    else
      v_i => v_c
    end if

    ! pass parent data to children .............................................

    ! parent extracts data and starts sending
    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), v = v_p)
      call send_buf(i) % Send_Start()
    end do

    ! child prepares buffer and starts receiving
    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), v = v_r)
      call recv_buf(i) % Recv_Start()
    end do

    ! parent finalizes sending
    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    ! child finalizes receiving and assigns data to
    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(v = v_r)
    end do

    ! interpolate parent data to child variable ................................

    select case (iop % mode)
    case(0)
      ! identity
      call SetArray(v_i, v_r)
    case(1)
      ! p-refinement
      call TPO_AAA(iop%A(:,:,1), v_r, v_i)
    case(2)
      ! h- or hp-refinement
      call TPO_1To8(iop%A, v_r, v_i)
    end select

    if (intermediate_data) then
      do k = 1, n_comp

        ! copy active element data
        call SetArray( v_c(:,:,:,1:child%n_elem_active,k) &
                     , v_i(:,:,:,1:child%n_elem_active,k) )

        ! extract frozen element data
        do e = child%n_elem_active + 1, child%n_elem
          i = 8*(child%element(e)%cluster_id - 1) + child%element(e)%cluster_oct
          v_c(:,:,:,e,k) = v_i(:,:,:,i,k)
        end do

      end do
    end if

    ! finalization .............................................................

    deallocate(send_buf, send_map)
    deallocate(recv_buf, recv_map)
    deallocate(v_r)

    if (intermediate_data) then
      deallocate(v_i)
    end if

  end subroutine ParentToChildInterpolation_X

  !=============================================================================

end module Parent_To_Child_Interpolation__3D
