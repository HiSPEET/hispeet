!> summary:   Divergence element operator based on handcrafted kernels (R)
!> author:    Erik Pfister, Joerg Stiller
!> date:      2020/03/20, revised 2021/06/09
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_R__Hand_RDP
  use Kind_Parameters, only: RWP => RDP
  use TPO__Div__3D_R__Gen
  implicit none
  private

  public :: TPO_Div_R_Hand

  interface TPO_Div_R_Hand
    module procedure TPO_Div_R_Hand_RWP
  end interface TPO_Div_R_Hand

contains

  !-----------------------------------------------------------------------------
  !> Divergence element operator based on handcrafted kernels (R)

  subroutine TPO_Div_R_Hand_RWP(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RWP), intent(in)  :: Ds(np,np)        !< standard diff matrix
    real(RWP), intent(in)  :: dx(3)            !< element extensions
    real(RWP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
    real(RWP), intent(out) :: v(np,np,np,ne)   !< element-wise divergence of u

    select case(np)
    case( 2)
      call TPO_Div_R_Hand__2  (ne, Ds, dx, u, v)
    case( 3)
      call TPO_Div_R_Hand__3  (ne, Ds, dx, u, v)
    case( 4)
      call TPO_Div_R_Hand__4  (ne, Ds, dx, u, v)
    case( 5)
      call TPO_Div_R_Hand__5  (ne, Ds, dx, u, v)
    case( 6)
      call TPO_Div_R_Hand__6  (ne, Ds, dx, u, v)
    case( 7)
      call TPO_Div_R_Hand__7  (ne, Ds, dx, u, v)
#ifndef __LIBXSMM__
    case( 8)
      call TPO_Div_R_Hand__8  (ne, Ds, dx, u, v)
    case( 9)
      call TPO_Div_R_Hand__9  (ne, Ds, dx, u, v)
    case(10)
      call TPO_Div_R_Hand__10 (ne, Ds, dx, u, v)
    case(11)
      call TPO_Div_R_Hand__11 (ne, Ds, dx, u, v)
    case(12)
      call TPO_Div_R_Hand__12 (ne, Ds, dx, u, v)
    case(13)
      call TPO_Div_R_Hand__13 (ne, Ds, dx, u, v)
    case(14)
      call TPO_Div_R_Hand__14 (ne, Ds, dx, u, v)
    case(15)
      call TPO_Div_R_Hand__15 (ne, Ds, dx, u, v)
    case(16)
      call TPO_Div_R_Hand__16 (ne, Ds, dx, u, v)
    case(17)
      call TPO_Div_R_Hand__17 (ne, Ds, dx, u, v)
    case(18)
      call TPO_Div_R_Hand__18 (ne, Ds, dx, u, v)
    case(19)
      call TPO_Div_R_Hand__19 (ne, Ds, dx, u, v)
    case(20)
      call TPO_Div_R_Hand__20 (ne, Ds, dx, u, v)
    case(21)
      call TPO_Div_R_Hand__21 (ne, Ds, dx, u, v)
    case(22)
      call TPO_Div_R_Hand__22 (ne, Ds, dx, u, v)
    case(23)
      call TPO_Div_R_Hand__23 (ne, Ds, dx, u, v)
    case(24)
      call TPO_Div_R_Hand__24 (ne, Ds, dx, u, v)
    case(25)
      call TPO_Div_R_Hand__25 (ne, Ds, dx, u, v)
    case(26)
      call TPO_Div_R_Hand__26 (ne, Ds, dx, u, v)
    case(27)
      call TPO_Div_R_Hand__27 (ne, Ds, dx, u, v)
    case(28)
      call TPO_Div_R_Hand__28 (ne, Ds, dx, u, v)
    case(29)
      call TPO_Div_R_Hand__29 (ne, Ds, dx, u, v)
    case(30)
      call TPO_Div_R_Hand__30 (ne, Ds, dx, u, v)
    case(31)
      call TPO_Div_R_Hand__31 (ne, Ds, dx, u, v)
    case(32)
      call TPO_Div_R_Hand__32 (ne, Ds, dx, u, v)
    case(33)
      call TPO_Div_R_Hand__33 (ne, Ds, dx, u, v)
#endif
    case default
      call TPO_Div_R_Gen(np, ne, Ds, dx, u, v)
    end select

  end subroutine TPO_Div_R_Hand_RWP

!===============================================================================
! operator routine

# ifdef __GFORTRAN__
#   define PASTE(base_name) base_name
#   define PROC(base_name,np) PASTE(base_name)np
# else
#  define PASTE(base_name,np) base_name ## np
#  define PROC(base_name,np) PASTE(base_name,np)
# endif

#define _NP_  2
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_

#ifndef __LIBXSMM__

#define _NP_  8
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__div__3d_r__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxQt.F"
#include "IxQtxI.F"
#include "QtxIxI.F"

!===============================================================================

end module TPO__Div__3D_R__Hand_RDP
