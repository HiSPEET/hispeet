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

!> summary:  Methods for transferring linked element face data of type real(RDP)
!> author:   Joerg Stiller
!> date:     2021/09/24
!===============================================================================

submodule(Element_Face_Transfer_Buffer__3D) SM_Real_RDP
  implicit none

   integer, parameter :: RK = RDP

contains

  !=============================================================================
  ! ElementFaceTransferBuffer_3D: constructor functions

  !-----------------------------------------------------------------------------
  !> New element face transfer buffer for a single real variable

  module function New_TransferBuffer_RDP_S(mesh, v) result(this)
    type(ElementFaceTransferBuffer_3D) :: this
    class(Mesh_3D), intent(in) :: mesh  !< mesh partition
    real(RDP), intent(in) :: v(:,:,:,:) !< element-face variable (np,np,6,*)

    this % np = size(v,1)
    this % nc = 1

    call Init_ElementFaceTransferBuffer_Shared_3D(this, mesh)

    if (this % nf <= 0) return ! no linked faces

    allocate(real(RDP) :: this % send % buf(this%np, this%np, this%nc, this%nf))
    allocate(real(RDP) :: this % recv % buf(this%np, this%np, this%nc, this%nf))

  end function New_TransferBuffer_RDP_S

  !-----------------------------------------------------------------------------
  !> New element face transfer buffer for an array of real variables

  module function New_TransferBuffer_RDP_A(mesh, v) result(this)
    type(ElementFaceTransferBuffer_3D) :: this
    class(Mesh_3D), intent(in) :: mesh    !< mesh partition
    real(RDP), intent(in) :: v(:,:,:,:,:) !< element-face variable (np,np,6,*,nc)

    this % np = size(v,1)
    this % nc = size(v,5)

    call Init_ElementFaceTransferBuffer_Shared_3D(this, mesh)

    if (this % nf <= 0) return ! no linked faces

    allocate(real(RDP) :: this % send % buf(this%np, this%np, this%nc, this%nf))
    allocate(real(RDP) :: this % recv % buf(this%np, this%np, this%nc, this%nf))

  end function New_TransferBuffer_RDP_A

  !=============================================================================
  ! ElementFaceTransferBuffer_3D: transfer procedures

  !----------------------------------------------------------------------------
  !> Send local data to ghosts and remote data -- real scalar

  module subroutine Transfer_RDP_S(this, mesh, v, tag)
    class(ElementFaceTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh       !< mesh partition
    real(RDP),      intent(in) :: v(:,:,:,:) !< element-face variable
    integer,        intent(in) :: tag        !< message tag

    if (this % nf <= 0) return ! no linked faces

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RDP_S

  !-----------------------------------------------------------------------------
  !> Send local data to ghosts and remote data -- real array

  module subroutine Transfer_RDP_A(this, mesh, v, tag)
    class(ElementFaceTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh         !< mesh partition
    real(RDP),      intent(in) :: v(:,:,:,:,:) !< element-face variable
    integer,        intent(in) :: tag          !< message tag

    if (this % nf <= 0) return ! no linked faces

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RDP_A

  !=============================================================================
  ! ElementFaceTransferBuffer_3D: merge procedures

  !---------------------------------------------------------------------------
  !> Complete receive and merge buffer into local data -- real(RDP) scalar

  module subroutine Merge_RDP_S(this, v)
    class(ElementFaceTransferBuffer_3D), intent(inout) :: this
    real(RDP), intent(inout) :: v(:,:,:,:) !< element-face variable

    if (this % nf <= 0) return ! no linked faces

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RDP_S'                      &
                , 'size(v,4) /= this%ne + this%ng'   &
                , 'Element_Face_Transfer_Buffer__3D' )
    end if
    !$omp end master

    call Merge_RX(this, v)

  end subroutine Merge_RDP_S

  !---------------------------------------------------------------------------
  !> Complete receive and merge buffer into local data -- real(RDP) array

  module subroutine Merge_RDP_A(this, v)
    class(ElementFaceTransferBuffer_3D), intent(inout) :: this
    real(RDP), intent(inout) :: v(:,:,:,:,:) !< element-face variable

    if (this % nf <= 0) return ! no linked faces

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RDP_A'                      &
                , 'size(v,4) /= this%ne + this%ng'   &
                , 'Element_Face_Transfer_Buffer__3D' )
    end if
    !$omp end master

    call Merge_RX(this, v)

  end subroutine Merge_RDP_A

  !=============================================================================
  ! Templated procedures Transfer_RX and Merge_RX

# include "element_face_transfer_buffer__3d__sm_real-tmpl.f"

  !=============================================================================

end submodule SM_Real_RDP
