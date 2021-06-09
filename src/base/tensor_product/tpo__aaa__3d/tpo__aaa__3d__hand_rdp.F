!> summary:   AxAxA kernel based on handcrafted kernels for real(RDP)
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__AAA__3D__Hand_RDP
  use Kind_Parameters, only: RWP => RDP
  use TPO__AAA__3D__Gen
  implicit none
  private

  public :: TPO_AAA_Hand

  interface TPO_AAA_Hand
    module procedure TPO_AAA_Hand_RWP
  end interface TPO_AAA_Hand

contains

  !-----------------------------------------------------------------------------
  !> Elliptic element operator based on handcrafted kernels (RCLI)

  subroutine TPO_AAA_Hand_RWP(A, u, v)
    real(RWP), intent(in)  :: A(:,:)     !< 1D operator
    real(RWP), intent(in)  :: u(:,:,:,:) !< operand
    real(RWP), intent(out) :: v(:,:,:,:) !< result

    integer :: na1, na2, ne

    na1 = size(A,1)
    na2 = size(A,2)
    ne  = size(u,4)

    ! 1:1 operators ............................................................

    if (na1 == na2) then

      select case(na1)
      case( 2)
        call TPO_AxAxA_Hand__2x2   (ne, A, u, v)
      case( 3)
        call TPO_AxAxA_Hand__3x3   (ne, A, u, v)
      case( 4)
        call TPO_AxAxA_Hand__4x4   (ne, A, u, v)
      case( 5)
        call TPO_AxAxA_Hand__5x5   (ne, A, u, v)
      case( 6)
        call TPO_AxAxA_Hand__6x6   (ne, A, u, v)
      case( 7)
        call TPO_AxAxA_Hand__7x7   (ne, A, u, v)
#ifndef __LIBXSMM__
      case( 8)
        call TPO_AxAxA_Hand__8x8   (ne, A, u, v)
      case( 9)
        call TPO_AxAxA_Hand__9x9   (ne, A, u, v)
      case(10)
        call TPO_AxAxA_Hand__10x10 (ne, A, u, v)
      case(11)
        call TPO_AxAxA_Hand__11x11 (ne, A, u, v)
      case(12)
        call TPO_AxAxA_Hand__12x12 (ne, A, u, v)
      case(13)
        call TPO_AxAxA_Hand__13x13 (ne, A, u, v)
      case(14)
        call TPO_AxAxA_Hand__14x14 (ne, A, u, v)
      case(15)
        call TPO_AxAxA_Hand__15x15 (ne, A, u, v)
      case(16)
        call TPO_AxAxA_Hand__16x16 (ne, A, u, v)
      case(17)
        call TPO_AxAxA_Hand__17x17 (ne, A, u, v)
      case(18)
        call TPO_AxAxA_Hand__18x18 (ne, A, u, v)
      case(19)
        call TPO_AxAxA_Hand__19x19 (ne, A, u, v)
      case(20)
        call TPO_AxAxA_Hand__20x20 (ne, A, u, v)
      case(21)
        call TPO_AxAxA_Hand__21x21 (ne, A, u, v)
      case(22)
        call TPO_AxAxA_Hand__22x22 (ne, A, u, v)
      case(23)
        call TPO_AxAxA_Hand__23x23 (ne, A, u, v)
      case(24)
        call TPO_AxAxA_Hand__24x24 (ne, A, u, v)
      case(25)
        call TPO_AxAxA_Hand__25x25 (ne, A, u, v)
      case(26)
        call TPO_AxAxA_Hand__26x26 (ne, A, u, v)
      case(27)
        call TPO_AxAxA_Hand__27x27 (ne, A, u, v)
      case(28)
        call TPO_AxAxA_Hand__28x28 (ne, A, u, v)
      case(29)
        call TPO_AxAxA_Hand__29x29 (ne, A, u, v)
      case(30)
        call TPO_AxAxA_Hand__30x30 (ne, A, u, v)
      case(31)
        call TPO_AxAxA_Hand__31x31 (ne, A, u, v)
      case(32)
        call TPO_AxAxA_Hand__32x32 (ne, A, u, v)
      case(33)
        call TPO_AxAxA_Hand__33x33 (ne, A, u, v)
#endif
      case default
        call TPO_AAA_Gen(A, u, v)
      end select

  ! 1:2 and 2:1 operators ......................................................

  else if (na1 == 2) then
    if (na2 == 3) then
      call TPO_AxAxA_Hand__2x3 (ne, A, u, v)
    else if (na2 == 4) then
      call TPO_AxAxA_Hand__2x4 (ne, A, u, v)
    else
      call TPO_AAA_Gen(A, u, v)
    end if

  else if (na1 == 3) then
    if (na2 == 4) then
      call TPO_AxAxA_Hand__3x4 (ne, A, u, v)
    else if (na2 == 5) then
      call TPO_AxAxA_Hand__3x5 (ne, A, u, v)
    else
      call TPO_AAA_Gen(A, u, v)
    end if

#ifndef __LIBXSMM__
  else if (na1 ==  4 .and. na2 ==  8) then
    call TPO_AxAxA_Hand__4x8 (ne, A, u, v)
  else if (na1 ==  5 .and. na2 ==  9) then
    call TPO_AxAxA_Hand__5x9 (ne, A, u, v)
  else if (na1 ==  8 .and. na2 == 16) then
    call TPO_AxAxA_Hand__8x16 (ne, A, u, v)
  else if (na1 ==  9 .and. na2 == 17) then
    call TPO_AxAxA_Hand__9x17 (ne, A, u, v)
  else if (na1 == 16 .and. na2 == 30) then
    call TPO_AxAxA_Hand__16x30 (ne, A, u, v)
  else if (na1 == 17 .and. na2 == 33) then
    call TPO_AxAxA_Hand__17x33 (ne, A, u, v)
#endif

  else if (na2 == 2) then
    if (na1 == 3) then
      call TPO_AxAxA_Hand__3x2 (ne, A, u, v)
    else if (na1 == 4) then
      call TPO_AxAxA_Hand__4x2 (ne, A, u, v)
    else
      call TPO_AAA_Gen(A, u, v)
    end if

  else if (na2 == 3) then
    if (na1 == 4) then
      call TPO_AxAxA_Hand__4x3 (ne, A, u, v)
    else if (na1 == 5) then
      call TPO_AxAxA_Hand__5x3 (ne, A, u, v)
    else
      call TPO_AAA_Gen(A, u, v)
    end if

#ifndef __LIBXSMM__
  else if (na1 ==  8 .and. na2 ==  4) then
    call TPO_AxAxA_Hand__8x4 (ne, A, u, v)
  else if (na1 ==  9 .and. na2 ==  5) then
    call TPO_AxAxA_Hand__9x5 (ne, A, u, v)
  else if (na1 == 16 .and. na2 ==  8) then
    call TPO_AxAxA_Hand__16x8 (ne, A, u, v)
  else if (na1 == 17 .and. na2 ==  9) then
    call TPO_AxAxA_Hand__17x9 (ne, A, u, v)
  else if (na1 == 30 .and. na2 == 16) then
    call TPO_AxAxA_Hand__30x16 (ne, A, u, v)
  else if (na1 == 33 .and. na2 == 17) then
    call TPO_AxAxA_Hand__33x17 (ne, A, u, v)
#endif

  else

    call TPO_AAA_Gen(A, u, v)

  end if

  end subroutine TPO_AAA_Hand_RWP


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
! 1:1 operators for na1 = na2 = 2:7

#define _NA1_  2
#define _NA2_  2
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  3
#define _NA2_  3
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  4
#define _NA2_  4
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  5
#define _NA2_  5
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  6
#define _NA2_  6
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  7
#define _NA2_  7
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

!-------------------------------------------------------------------------------
! 1:1 operators for na1 = na2 = 8:33

#ifndef __LIBXSMM__

#define _NA1_  8
#define _NA2_  8
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  9
#define _NA2_  9
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  10
#define _NA2_  10
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  11
#define _NA2_  11
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  12
#define _NA2_  12
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  13
#define _NA2_  13
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  14
#define _NA2_  14
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  15
#define _NA2_  15
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  16
#define _NA2_  16
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  17
#define _NA2_  17
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  18
#define _NA2_  18
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  19
#define _NA2_  19
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  20
#define _NA2_  20
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  21
#define _NA2_  21
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  22
#define _NA2_  22
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  23
#define _NA2_  23
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  24
#define _NA2_  24
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  25
#define _NA2_  25
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  26
#define _NA2_  26
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  27
#define _NA2_  27
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  28
#define _NA2_  28
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  29
#define _NA2_  29
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  30
#define _NA2_  30
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  31
#define _NA2_  31
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  32
#define _NA2_  32
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  33
#define _NA2_  33
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

#endif

!-------------------------------------------------------------------------------
! 1:2 and 2:1 operators up to max(na1,na2) = 7

#define _NA1_  2
#define _NA2_  3
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  2
#define _NA2_  4
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  3
#define _NA2_  4
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  3
#define _NA2_  5
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  3
#define _NA2_  2
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  4
#define _NA2_  2
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  4
#define _NA2_  3
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  5
#define _NA2_  3
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

!-------------------------------------------------------------------------------
! remaining 1:2 and 2:1 operators up to max(na1,na2) = 33

#ifndef __LIBXSMM__

#define _NA1_  4
#define _NA2_  8
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  5
#define _NA2_  9
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  8
#define _NA2_ 16
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  9
#define _NA2_ 17
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 16
#define _NA2_ 30
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 17
#define _NA2_ 33
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  8
#define _NA2_  4
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_  9
#define _NA2_  5
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 16
#define _NA2_  8
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 17
#define _NA2_  9
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 30
#define _NA2_ 16
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_
#define _NA1_ 33
#define _NA2_ 17
#include "tpo__aaa__3d__hand_kernel.f"
#undef _NA1_
#undef _NA2_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxAt.F"
#include "IxBtxI.F"
#include "CtxIxI.F"

!===============================================================================

end module TPO__AAA__3D__Hand_RDP
