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

submodule(CART__TPO_Diffusion_CI) MP_AssignProcEVec
  use Kind_Parameters,  only: RNP
  implicit none

  ! these parameters will be moved into a separate module
  integer, parameter :: LVEC    = 4
  integer, parameter :: LVEC_M1 = LVEC - 1

contains

!===============================================================================
! Selector

module subroutine AssignProcEVec(np, Proc)
  integer, intent(in) :: np  !< number of points per direction
  procedure(TPO_Diffusion_CI_Proc), pointer, intent(out) :: Proc

  select case(np)
  case(2)
    Proc => CART__TPO_Diffusion_CI__EVec_2
  ! ... ausfüllen ...
  case(15)
    Proc => CART__TPO_Diffusion_CI__EVec_15
  case default
    Proc => null()
  end select

end subroutine AssignProcEVec

!===============================================================================
! Parametrized routines

#define _PROC(NP)   CART__TPO_Diffusion_CI__EVec_ ## NP
#define _PROC_(NP)  _PROC(NP)

#define  _NP_  2
#include "cart__TPO_Diffusion_ci__evec_par.ft"
#undef   _NP_

! ... ausfüllen ... bis 15 sollte erstmal reichen

#define  _NP_  15
#include "cart__TPO_Diffusion_ci__evec_par.ft"
#undef   _NP_

!===============================================================================

end submodule MP_AssignProcEVec
