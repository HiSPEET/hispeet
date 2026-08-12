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

!> summary:  Initialization of HDF5 and provison of custom intrinsic types
!> author:   Joerg Stiller
!> date:     2024/06/04
!==============================================================================

module HDF5_Binding
  use Kind_Parameters, only: IXS, IXL, RSP, RDP, RHP, RNP
  use HDF5
  use H5LT
  public

  !-----------------------------------------------------------------------------
  ! HDF5 data types, naming follows HDF5 conventions

  integer(HID_T), save :: &
    H5T_INTEGER,     & !< HDF5 data type matching integer
    H5T_INTEGER_IXS, & !< HDF5 data type matching integer(IXS)
    H5T_INTEGER_IXL, & !< HDF5 data type matching integer(IXL)
    H5T_REAL_RSP,    & !< HDF5 data type matching real(RDP)
    H5T_REAL_RDP,    & !< HDF5 data type matching real(RDP)
    H5T_REAL_RHP,    & !< HDF5 data type matching real(RHP)
    H5T_REAL_RNP,    & !< HDF5 data type matching real(RNP)
    H5T_CHARACTER,   & !< HDF5 data type matching character
    H5T_LOGICAL        !< HDF5 data type matching logical

  !-----------------------------------------------------------------------------
  ! Private variables

  logical, private :: initialized = .false.

contains

  !-----------------------------------------------------------------------------
  !> Initializes the HDF5 data types defined in the module
  !>
  !> HDF5 does not support logical yet. H5T_INTEGER is used as a workaround.
  !> This works with GCC and Intel, but may fail with other compilers.

  subroutine Init_HDF5_Binding()

    integer :: err

    if (initialized) return

    ! safeguard, there are no damaging side effects in calling it more than once
    call H5open_f(err)

    H5T_INTEGER     = H5T_NATIVE_INTEGER
    H5T_INTEGER_IXS = H5kind_to_type(IXS, H5_INTEGER_KIND)
    H5T_INTEGER_IXL = H5kind_to_type(IXL, H5_INTEGER_KIND)

    H5T_REAL_RSP    = H5T_NATIVE_REAL
    H5T_REAL_RDP    = H5T_NATIVE_DOUBLE
    H5T_REAL_RHP    = H5kind_to_type(RHP, H5_REAL_KIND)
    H5T_REAL_RNP    = H5kind_to_type(RNP, H5_REAL_KIND)

    H5T_CHARACTER   = H5T_NATIVE_CHARACTER
    H5T_LOGICAL     = H5T_INTEGER

    initialized = .true.

  end subroutine Init_HDF5_Binding

  !============================================================================

end module HDF5_Binding
