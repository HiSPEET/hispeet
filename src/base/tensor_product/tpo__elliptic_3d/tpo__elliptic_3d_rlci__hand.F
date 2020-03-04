!> summary:   Elliptic element operator based on handcrafted kernels (RCLI)
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic_3d_RLCI__Hand
  use Kind_Parameters, only: RNP
  use TPO__Elliptic_3d_RLCI__Gen
  implicit none
  private

  public :: TPO_Elliptic_RLCI_Hand

contains

  !-----------------------------------------------------------------------------
  !> Elliptic element operator based on handcrafted kernels (RCLI)

  subroutine TPO_Elliptic_RLCI_Hand(np, ne, Ms, Ls, lambda, nu, dx, u, v)
    integer,   intent(in)  :: np             !< num points per direction
    integer,   intent(in)  :: ne             !< num elements
    real(RNP), intent(in)  :: Ms(np)         !< standard mass matrix
    real(RNP), intent(in)  :: Ls(np,np)      !< standard stiffness matrix
    real(RNP), intent(in)  :: lambda         !< Helmholtz parameter λ
    real(RNP), intent(in)  :: nu             !< diffusivity nu
    real(RNP), intent(in)  :: dx(3)          !< element extensions
    real(RNP), intent(in)  :: u(np,np,np,ne) !< operand
    real(RNP), intent(out) :: v(np,np,np,ne) !< result

    select case(np)
    case( 2)
      call TPO_Elliptic_RLCI_Hand__2  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 3)
      call TPO_Elliptic_RLCI_Hand__3  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 4)
      call TPO_Elliptic_RLCI_Hand__4  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 5)
      call TPO_Elliptic_RLCI_Hand__5  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 6)
      call TPO_Elliptic_RLCI_Hand__6  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 7)
      call TPO_Elliptic_RLCI_Hand__7  (ne, Ms, Ls, lambda, nu, dx, u, v)
#ifndef LIBXSMM
    case( 8)
      call TPO_Elliptic_RLCI_Hand__8  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case( 9)
      call TPO_Elliptic_RLCI_Hand__9  (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(10)
      call TPO_Elliptic_RLCI_Hand__10 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(11)
      call TPO_Elliptic_RLCI_Hand__11 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(12)
      call TPO_Elliptic_RLCI_Hand__12 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(13)
      call TPO_Elliptic_RLCI_Hand__13 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(14)
      call TPO_Elliptic_RLCI_Hand__14 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(15)
      call TPO_Elliptic_RLCI_Hand__15 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(16)
      call TPO_Elliptic_RLCI_Hand__16 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(17)
      call TPO_Elliptic_RLCI_Hand__17 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(18)
      call TPO_Elliptic_RLCI_Hand__18 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(19)
      call TPO_Elliptic_RLCI_Hand__19 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(20)
      call TPO_Elliptic_RLCI_Hand__20 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(21)
      call TPO_Elliptic_RLCI_Hand__21 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(22)
      call TPO_Elliptic_RLCI_Hand__22 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(23)
      call TPO_Elliptic_RLCI_Hand__23 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(24)
      call TPO_Elliptic_RLCI_Hand__24 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(25)
      call TPO_Elliptic_RLCI_Hand__25 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(26)
      call TPO_Elliptic_RLCI_Hand__26 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(27)
      call TPO_Elliptic_RLCI_Hand__27 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(28)
      call TPO_Elliptic_RLCI_Hand__28 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(29)
      call TPO_Elliptic_RLCI_Hand__29 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(30)
      call TPO_Elliptic_RLCI_Hand__30 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(31)
      call TPO_Elliptic_RLCI_Hand__31 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(32)
      call TPO_Elliptic_RLCI_Hand__32 (ne, Ms, Ls, lambda, nu, dx, u, v)
    case(33)
      call TPO_Elliptic_RLCI_Hand__33 (ne, Ms, Ls, lambda, nu, dx, u, v)
#endif
    case default
      call TPO_Elliptic_RLCI_Gen(np, ne, Ms, Ls, lambda, nu, dx, u, v)
    end select

  end subroutine TPO_Elliptic_RLCI_Hand

  !-----------------------------------------------------------------------------
  !> Projection to element mass matrix and computation of the Helmholtz term
  !>
  !>   * k:  joined with i,j and flattened  (f)
  !>   * j:  joined with i,k and flattened  (f)
  !>   * i:  joined with j,k and flattened  (f)

  subroutine SetOperands(np_e, lambda, M, u, M_u, v)
    !$acc routine vector
    integer,   intent(in)  :: np_e
    real(RNP), intent(in)  :: lambda
    real(RNP), intent(in)  :: M   (np_e)
    real(RNP), intent(in)  :: u   (np_e)
    real(RNP), intent(out) :: M_u (np_e)
    real(RNP), intent(out) :: v   (np_e)

    integer :: l

    !$acc loop vector
    !DIR$ SIMD
    do l = 1, np_e
      M_u(l) = M(l) * u(l)
      v(l) = lambda * M_u(l)
    end do

  end subroutine SetOperands

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
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_

#ifndef LIBXSMM

#define _NP_  8
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__elliptic_3d_rlci__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "../primitives/IxIxQt/IxIxQt.F"
#include "../primitives/IxQtxI/IxQtxI.F"
#include "../primitives/QtxIxI/QtxIxI.F"

!===============================================================================

end module TPO__Elliptic_3d_RLCI__Hand
