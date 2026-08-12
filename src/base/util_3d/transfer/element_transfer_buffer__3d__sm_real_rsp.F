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

!> summary:  Methods for transferring linked element data of type real(RSP)
!> author:   Joerg Stiller
!> date:     2021/02/20
!===============================================================================

submodule(Element_Transfer_Buffer__3D) SM_Real_RSP
  implicit none

  integer, parameter :: RK = RSP

contains

  !=============================================================================
  ! ElementTransferBuffer_3D: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single real variable

  module function New_TransferBuffer_RSP_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    class(Mesh_3D),    intent(in) :: mesh       !< mesh partition
    real(RSP),         intent(in) :: v(:,:,:,:) !< mesh variable
    integer, optional, intent(in) :: nl(3)      !< number of layers

    integer :: np(3)

    if (mesh % part   < 0) return
    if (mesh % n_link < 1) return

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = 1

    allocate(real(RSP) :: this % master % buf(size(this % master % node)))
    allocate(real(RSP) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_RSP_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of real variables

  module function New_TransferBuffer_RSP_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    class(Mesh_3D),    intent(in) :: mesh         !< mesh partition
    real(RSP),         intent(in) :: v(:,:,:,:,:) !< mesh variables
    integer, optional, intent(in) :: nl(3)        !< number of layers

    integer :: np(3), nc

    if (mesh % part   < 0) return
    if (mesh % n_link < 1) return

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = nc

    allocate(real(RSP) :: this % master % buf(size(this % master % node) * nc))
    allocate(real(RSP) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_RSP_A

  !=============================================================================
  ! ElementTransferBuffer_3D: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real scalar

  module subroutine Transfer_RSP_S(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh       !< mesh partition
    real(RSP),      intent(in) :: v(:,:,:,:) !< mesh variable
    integer,        intent(in) :: tag        !< message tag

    if (this % ng < 1) return

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RSP_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real array

  module subroutine Transfer_RSP_A(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh         !< mesh partition
    real(RSP),      intent(in) :: v(:,:,:,:,:) !< mesh variable
    integer,        intent(in) :: tag          !< message tag

    if (this % ng < 1) return

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RSP_A

  !=============================================================================
  ! ElementTransferBuffer_3D: merge procedures

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real scalar
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng)`, where
  !>
  !>    *  `np` is the number of points per direction, as in this%np
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements

  module subroutine Merge_RSP_S(this, v)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    real(RSP), intent(inout) :: v(:,:,:,:)  !< mesh variable

    if (this % ng < 1) return

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RSP_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_RX(this, v)

  end subroutine Merge_RSP_S

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real array
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng, nc)`, where
  !>
  !>    *  `np` is the number of points per direction, as in `this%np`
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements
  !>    *  `nc` is the number of components, as in `this%nc`

  module subroutine Merge_RSP_A(this, v)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    real(RSP), intent(inout) :: v(:,:,:,:,:)  !< mesh variable

    if (this % ng < 1) return

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RSP_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_RX(this, v)

  end subroutine Merge_RSP_A

  !=============================================================================
  ! Templated procedures Transfer_RX and Merge_RX

# include "element_transfer_buffer__3d__sm_real-tmpl.f"

  !=============================================================================

end submodule SM_Real_RSP
