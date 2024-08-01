!> summary:  3d regular 1:8 interpolation with handcrafted kernels for real(RDP)
!> author:   Joerg Stiller
!> date:     2024/07/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__1To8__3D__Hand_RDP
  use Kind_Parameters, only: RWP => RDP
  use TPO__1To8__3D__Gen
  implicit none
  private

  public :: TPO_1To8_Hand

  interface TPO_1To8_Hand
    module procedure TPO_1To8_Hand_RWP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> element diffusion operator based on handcrafted kernels (RCLI)

  subroutine TPO_1To8_Hand_RWP(A, u, v)
    real(RWP), intent(in)  :: A(:,:,:)   !< 1D operator
    real(RWP), intent(in)  :: u(:,:,:,:) !< operand
    real(RWP), intent(out) :: v(:,:,:,:) !< result

    integer :: na1, na2, ne

    na1 = size(A,1)
    na2 = size(A,2)
    ne  = size(u,4)

    if (na1 == na2) then

      select case(na1)
      case( 2)
        call TPO_1To8_Hand__2x2   (ne, A, u, v)
      case( 3)
        call TPO_1To8_Hand__3x3   (ne, A, u, v)
      case( 4)
        call TPO_1To8_Hand__4x4   (ne, A, u, v)
      case( 5)
        call TPO_1To8_Hand__5x5   (ne, A, u, v)
      case( 6)
        call TPO_1To8_Hand__6x6   (ne, A, u, v)
      case( 7)
        call TPO_1To8_Hand__7x7   (ne, A, u, v)
#ifndef __LIBXSMM__
      case( 8)
        call TPO_1To8_Hand__8x8   (ne, A, u, v)
      case( 9)
        call TPO_1To8_Hand__9x9   (ne, A, u, v)
      case(10)
        call TPO_1To8_Hand__10x10 (ne, A, u, v)
      case(11)
        call TPO_1To8_Hand__11x11 (ne, A, u, v)
      case(12)
        call TPO_1To8_Hand__12x12 (ne, A, u, v)
      case(13)
        call TPO_1To8_Hand__13x13 (ne, A, u, v)
      case(14)
        call TPO_1To8_Hand__14x14 (ne, A, u, v)
      case(15)
        call TPO_1To8_Hand__15x15 (ne, A, u, v)
      case(16)
        call TPO_1To8_Hand__16x16 (ne, A, u, v)
      case(17)
        call TPO_1To8_Hand__17x17 (ne, A, u, v)
      case(18)
        call TPO_1To8_Hand__18x18 (ne, A, u, v)
      case(19)
        call TPO_1To8_Hand__19x19 (ne, A, u, v)
      case(20)
        call TPO_1To8_Hand__20x20 (ne, A, u, v)
      case(21)
        call TPO_1To8_Hand__21x21 (ne, A, u, v)
      case(22)
        call TPO_1To8_Hand__22x22 (ne, A, u, v)
      case(23)
        call TPO_1To8_Hand__23x23 (ne, A, u, v)
      case(24)
        call TPO_1To8_Hand__24x24 (ne, A, u, v)
      case(25)
        call TPO_1To8_Hand__25x25 (ne, A, u, v)
      case(26)
        call TPO_1To8_Hand__26x26 (ne, A, u, v)
      case(27)
        call TPO_1To8_Hand__27x27 (ne, A, u, v)
      case(28)
        call TPO_1To8_Hand__28x28 (ne, A, u, v)
      case(29)
        call TPO_1To8_Hand__29x29 (ne, A, u, v)
      case(30)
        call TPO_1To8_Hand__30x30 (ne, A, u, v)
      case(31)
        call TPO_1To8_Hand__31x31 (ne, A, u, v)
      case(32)
        call TPO_1To8_Hand__32x32 (ne, A, u, v)
      case(33)
        call TPO_1To8_Hand__33x33 (ne, A, u, v)
#endif
      case default
        call TPO_1To8_Gen(A, u, v)
      end select

  else

    call TPO_1To8_Gen(A, u, v)

  end if

  end subroutine TPO_1To8_Hand_RWP


!===============================================================================
! operator routines

# ifdef __GFORTRAN__
#   define PASTE(base_name) base_name
#   define PROC(base_name,na1,na2) PASTE(base_name)PASTE(na1)PASTE(x)na2
# else
#  define PASTE(base_name,na1,na2) base_name ## na1 ## x ##na2
#  define PROC(base_name,na1,na2) PASTE(base_name,na1,na2)
# endif

!-------------------------------------------------------------------------------
! operators for na1 = na2 = 2:7

#define _NA1_  2
#define _NA2_  2
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  3
#define _NA2_  3
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  4
#define _NA2_  4
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  5
#define _NA2_  5
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  6
#define _NA2_  6
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  7
#define _NA2_  7
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

!-------------------------------------------------------------------------------
! 1:1 operators for na1 = na2 = 8:33

#ifndef __LIBXSMM__

#define _NA1_  8
#define _NA2_  8
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  9
#define _NA2_  9
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  10
#define _NA2_  10
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  11
#define _NA2_  11
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  12
#define _NA2_  12
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  13
#define _NA2_  13
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  14
#define _NA2_  14
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  15
#define _NA2_  15
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  16
#define _NA2_  16
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  17
#define _NA2_  17
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  18
#define _NA2_  18
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  19
#define _NA2_  19
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  20
#define _NA2_  20
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  21
#define _NA2_  21
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  22
#define _NA2_  22
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  23
#define _NA2_  23
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  24
#define _NA2_  24
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  25
#define _NA2_  25
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  26
#define _NA2_  26
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  27
#define _NA2_  27
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  28
#define _NA2_  28
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  29
#define _NA2_  29
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  30
#define _NA2_  30
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  31
#define _NA2_  31
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  32
#define _NA2_  32
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  33
#define _NA2_  33
#include "tpo__1To8__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxAt.F"
#include "IxBtxI.F"
#include "CtxIxI.F"

!===============================================================================

end module TPO__1To8__3D__Hand_RDP
