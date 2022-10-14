!> summary:   Gradient element operator based on handcrafted kernels (D)
!> author:    Jerome Michel
!> date:      2021/10/12
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D_D__Hand_RDP
  use Kind_Parameters, only: RWP => RDP
  use TPO__Grad__3D_D__Gen
  implicit none
  private

  public :: TPO_Grad_D_Hand

  interface TPO_Grad_D_Hand
    module procedure TPO_Grad_D_Hand_RWP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Gradient element operator based on handcrafted kernels (D)

  subroutine TPO_Grad_D_Hand_RWP(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    integer,   intent(in)  :: np                  !< num points per direction
    integer,   intent(in)  :: ne                  !< num elements
    real(RWP), intent(in)  :: Ms(np)              !< 1D standard mass matrix
    real(RWP), intent(in)  :: Ds(np,np)           !< 1D standard diff matrix
    real(RWP), intent(in)  :: Jd(np,np,np,ne)     !< element Jacobian determinant
    real(RWP), intent(in)  :: Ji(np,np,np,ne,3,3) !< inverse Jacobi matrix
    real(RWP), intent(in)  :: a(np,np,6,ne)       !< face area coefficients
    real(RWP), intent(in)  :: n(np,np,6,ne,3)     !< face unit normal vectors
    real(RWP), intent(in)  :: u(np,np,np,ne)      !< 3D scalar field u
    real(RWP), intent(in)  :: up(np,np,6,ne)      !< exterior traces u⁺ at faces
    real(RWP), intent(out) :: v(np,np,np,ne,3)    !< element-wise gradient of u

    optional :: Ms, Jd, a, n, up

    select case(np)
    case( 2)
      call TPO_Grad_D_Hand__2  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 3)
      call TPO_Grad_D_Hand__3  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 4)
      call TPO_Grad_D_Hand__4  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 5)
      call TPO_Grad_D_Hand__5  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 6)
      call TPO_Grad_D_Hand__6  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 7)
      call TPO_Grad_D_Hand__7  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
#ifndef __LIBXSMM__
    case( 8)
      call TPO_Grad_D_Hand__8  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case( 9)
      call TPO_Grad_D_Hand__9  (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(10)
      call TPO_Grad_D_Hand__10 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(11)
      call TPO_Grad_D_Hand__11 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(12)
      call TPO_Grad_D_Hand__12 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(13)
      call TPO_Grad_D_Hand__13 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(14)
      call TPO_Grad_D_Hand__14 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(15)
      call TPO_Grad_D_Hand__15 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(16)
      call TPO_Grad_D_Hand__16 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(17)
      call TPO_Grad_D_Hand__17 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(18)
      call TPO_Grad_D_Hand__18 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(19)
      call TPO_Grad_D_Hand__19 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(20)
      call TPO_Grad_D_Hand__20 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(21)
      call TPO_Grad_D_Hand__21 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(22)
      call TPO_Grad_D_Hand__22 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(23)
      call TPO_Grad_D_Hand__23 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(24)
      call TPO_Grad_D_Hand__24 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(25)
      call TPO_Grad_D_Hand__25 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(26)
      call TPO_Grad_D_Hand__26 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(27)
      call TPO_Grad_D_Hand__27 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(28)
      call TPO_Grad_D_Hand__28 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(29)
      call TPO_Grad_D_Hand__29 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(30)
      call TPO_Grad_D_Hand__30 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(31)
      call TPO_Grad_D_Hand__31 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(32)
      call TPO_Grad_D_Hand__32 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    case(33)
      call TPO_Grad_D_Hand__33 (ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
#endif
    case default
      call TPO_Grad_D_Gen(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    end select

  end subroutine TPO_Grad_D_Hand_RWP

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
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_

#ifndef __LIBXSMM__

#define _NP_  8
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__grad__3d_d__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxQt.F"
#include "IxQtxI.F"
#include "QtxIxI.F"

!===============================================================================

end module TPO__Grad__3D_D__Hand_RDP
