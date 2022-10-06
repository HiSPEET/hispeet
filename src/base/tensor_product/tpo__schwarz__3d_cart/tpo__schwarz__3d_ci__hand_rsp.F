!> summary:   3d isotropic Schwarz operator based on handcrafted kernels (CI)
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/06/02
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_CI__Hand_RSP
  use Kind_Parameters, only: RWP => RSP
  use TPO__Schwarz__3D_CI__Gen
  implicit none
  private

  public :: TPO_Schwarz_CI_Hand

  interface TPO_Schwarz_CI_Hand
    module procedure TPO_Schwarz_CI_Hand_RWP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Isotropic schwarz operator based on handcrafted kernels (I)

  subroutine TPO_Schwarz_CI_Hand_RWP(S, W, cfg, D_inv, f, u)
    real(RWP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RWP), intent(in)  :: W(:,:)         !< 1D weights
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs
    real(RWP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvals
    real(RWP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RWP), intent(out) :: u(:,:,:,:)     !< solution

    integer :: nc, nd, np

    np = size(S,1)
    nc = size(S,3)
    nd = size(f,4)

    select case(np)
    case( 2)
      call TPO_Schwarz_CI_Hand__2  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 3)
      call TPO_Schwarz_CI_Hand__3  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 4)
      call TPO_Schwarz_CI_Hand__4  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 5)
      call TPO_Schwarz_CI_Hand__5  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 6)
      call TPO_Schwarz_CI_Hand__6  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 7)
      call TPO_Schwarz_CI_Hand__7  (nc, nd, S, W, cfg, D_inv, f, u)
#ifndef __LIBXSMM__
    case( 8)
      call TPO_Schwarz_CI_Hand__8  (nc, nd, S, W, cfg, D_inv, f, u)
    case( 9)
      call TPO_Schwarz_CI_Hand__9  (nc, nd, S, W, cfg, D_inv, f, u)
    case(10)
      call TPO_Schwarz_CI_Hand__10 (nc, nd, S, W, cfg, D_inv, f, u)
    case(11)
      call TPO_Schwarz_CI_Hand__11 (nc, nd, S, W, cfg, D_inv, f, u)
    case(12)
      call TPO_Schwarz_CI_Hand__12 (nc, nd, S, W, cfg, D_inv, f, u)
    case(13)
      call TPO_Schwarz_CI_Hand__13 (nc, nd, S, W, cfg, D_inv, f, u)
    case(14)
      call TPO_Schwarz_CI_Hand__14 (nc, nd, S, W, cfg, D_inv, f, u)
    case(15)
      call TPO_Schwarz_CI_Hand__15 (nc, nd, S, W, cfg, D_inv, f, u)
    case(16)
      call TPO_Schwarz_CI_Hand__16 (nc, nd, S, W, cfg, D_inv, f, u)
    case(17)
      call TPO_Schwarz_CI_Hand__17 (nc, nd, S, W, cfg, D_inv, f, u)
    case(18)
      call TPO_Schwarz_CI_Hand__18 (nc, nd, S, W, cfg, D_inv, f, u)
    case(19)
      call TPO_Schwarz_CI_Hand__19 (nc, nd, S, W, cfg, D_inv, f, u)
    case(20)
      call TPO_Schwarz_CI_Hand__20 (nc, nd, S, W, cfg, D_inv, f, u)
    case(21)
      call TPO_Schwarz_CI_Hand__21 (nc, nd, S, W, cfg, D_inv, f, u)
    case(22)
      call TPO_Schwarz_CI_Hand__22 (nc, nd, S, W, cfg, D_inv, f, u)
    case(23)
      call TPO_Schwarz_CI_Hand__23 (nc, nd, S, W, cfg, D_inv, f, u)
    case(24)
      call TPO_Schwarz_CI_Hand__24 (nc, nd, S, W, cfg, D_inv, f, u)
    case(25)
      call TPO_Schwarz_CI_Hand__25 (nc, nd, S, W, cfg, D_inv, f, u)
    case(26)
      call TPO_Schwarz_CI_Hand__26 (nc, nd, S, W, cfg, D_inv, f, u)
    case(27)
      call TPO_Schwarz_CI_Hand__27 (nc, nd, S, W, cfg, D_inv, f, u)
    case(28)
      call TPO_Schwarz_CI_Hand__28 (nc, nd, S, W, cfg, D_inv, f, u)
    case(29)
      call TPO_Schwarz_CI_Hand__29 (nc, nd, S, W, cfg, D_inv, f, u)
    case(30)
      call TPO_Schwarz_CI_Hand__30 (nc, nd, S, W, cfg, D_inv, f, u)
    case(31)
      call TPO_Schwarz_CI_Hand__31 (nc, nd, S, W, cfg, D_inv, f, u)
    case(32)
      call TPO_Schwarz_CI_Hand__32 (nc, nd, S, W, cfg, D_inv, f, u)
    case(33)
      call TPO_Schwarz_CI_Hand__33 (nc, nd, S, W, cfg, D_inv, f, u)
#endif
    case default
      call TPO_Schwarz_CI_Gen(S, W, cfg, D_inv, f, u)
    end select

  end subroutine TPO_Schwarz_CI_Hand_RWP

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
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_

#ifndef __LIBXSMM__

#define _NP_  8
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__schwarz__3d_ci__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxQt.F"
#include "IxQtxI.F"
#include "QtxIxI.F"

!===============================================================================

end module TPO__Schwarz__3D_CI__Hand_RSP
