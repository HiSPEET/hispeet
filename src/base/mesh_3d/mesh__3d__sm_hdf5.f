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

!> summary:  HDF5 binding to mesh attributes
!> author:   Joerg Stiller, Erik Pfister, Moritz Kreuseler
!> date:     2024/06/05
!===============================================================================

submodule(Mesh__3D) SM_HDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

  !-----------------------------------------------------------------------------
  !> HDF5 datatype for static components of MeshAttributes_3D

  integer(HID_T) :: H5T_MeshAttributes = -1

contains

  !-----------------------------------------------------------------------------
  !> Get HDF5 datatype for essential static components of MeshAttributes_3D

  module subroutine Get_H5T_MeshAttributes_3D( H5T_MeshAttributes_3D )
    integer(HID_T), intent(out) :: H5T_MeshAttributes_3D

    !$omp master

    if (H5T_MeshAttributes < 0) then
      call Init_H5T_MeshAttributes()
    end if

    H5T_MeshAttributes_3D = H5T_MeshAttributes

    !$omp end master

  end subroutine Get_H5T_MeshAttributes_3D

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_MeshAttributes

  subroutine Init_H5T_MeshAttributes()

    type(MeshAttributes_3D), target :: attributes(2)
    integer(SIZE_T) :: offset
    integer(HID_T)  :: tid
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_MeshAttributes, err)

    ! insert n_bound
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%n_bound))
    call H5Tinsert_f(H5T_MeshAttributes, 'n_bound', offset, H5T_INTEGER, err)

    ! insert n_parts
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%n_parts))
    call H5Tinsert_f(H5T_MeshAttributes, 'n_parts', offset, H5T_INTEGER, err)

    ! insert p_geom
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%p_geom))
    call H5Tinsert_f(H5T_MeshAttributes, 'p_geom', offset, H5T_INTEGER, err)

    ! insert structured
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%structured))
    call H5Tinsert_f(H5T_MeshAttributes, 'structured', offset, H5T_LOGICAL, err)

    ! insert regular
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%regular))
    call H5Tinsert_f(H5T_MeshAttributes, 'regular', offset, H5T_LOGICAL, err)

    ! insert dx
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%dx(1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 1, int([3], HSIZE_T), tid, err)
    call H5Tinsert_f(H5T_MeshAttributes, 'dx', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert is_root
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%is_root))
    call H5Tinsert_f(H5T_MeshAttributes, 'is_root', offset, H5T_LOGICAL, err)

    ! insert is_top
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%is_top))
    call H5Tinsert_f(H5T_MeshAttributes, 'is_top', offset, H5T_LOGICAL, err)

    ! insert has_sfc
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%has_sfc))
    call H5Tinsert_f(H5T_MeshAttributes, 'has_sfc', offset, H5T_LOGICAL, err)

    ! insert refinement
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%refinement))
    call H5Tinsert_f(H5T_MeshAttributes, 'refinement', offset, H5T_CHARACTER, err)

  end subroutine Init_H5T_MeshAttributes

  !=============================================================================

end submodule SM_HDF5
