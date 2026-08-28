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

!> summary:  Multilevel SDC variable for Dahlquist
!> author:   Erik Pfister, Joerg Stiller
!> date:     2023/10/08
!===============================================================================

module DQ__MLSDC__Variable__1D
  use Kind_Parameters
  use DQ__MLSDC__1D

  implicit none
  private

  public :: DQ_MLSDC_Variable_1D

  !-----------------------------------------------------------------------------
  !> 1D multilevel space-time variable

  type DQ_MLSDC_Variable_1D
    type(LevelVariable), allocatable :: level(:)
  end type DQ_MLSDC_Variable_1D

  ! constructor
  interface DQ_MLSDC_Variable_1D
    procedure New_DQ_MLSDC_Variable_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Structure representing one level of a multilevel space-time variable

  type LevelVariable
    complex(RNP), allocatable :: val(:,:) !< variable values
  end type LevelVariable

contains

  !-----------------------------------------------------------------------------
  !> Construction of a new MLSDC variable
  !>
  !> The dimensions and, if not specified, the number of components are adopted
  !> from the given MLSDC data structure.

  function New_DQ_MLSDC_Variable_1D(mlsdc) result(var)
    class(DQ_MLSDC_1D), intent(in) :: mlsdc !< MLSDC data structure
    type(DQ_MLSDC_Variable_1D) :: var

    integer :: l, nl

    nl = size(mlsdc % level)

    allocate(var % level(nl))

    do l = 1, nl
      associate( mt => mlsdc % level(l) % m_time  &
               , nt => mlsdc % level(l) % n_time  )
        allocate( var % level(l) % val(0:mt,1:nt) &
                , source = cmplx(0.0, 0.0, kind=RNP) )
      end associate
    end do

  end function New_DQ_MLSDC_Variable_1D

  !=============================================================================

end module DQ__MLSDC__Variable__1D
