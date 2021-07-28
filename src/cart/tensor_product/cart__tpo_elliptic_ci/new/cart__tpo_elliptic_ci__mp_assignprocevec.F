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
