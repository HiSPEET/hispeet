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

!> summary:   3d isotropic Schwarz operator based on handcrafted kernels (CI)
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/06/02
!===============================================================================

module TPO__Schwarz__3D_I__Hand_RSP
  use Kind_Parameters, only: RWP => RSP
  use TPO__Schwarz__3D_I__Gen
  implicit none
  private

  public :: TPO_Schwarz_I_Hand

  interface TPO_Schwarz_I_Hand
    module procedure TPO_Schwarz_I_Hand_RWP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Isotropic schwarz operator based on handcrafted kernels (I)

  subroutine TPO_Schwarz_I_Hand_RWP(S, V, W, g, cfg, lambda, nu, f, u)
    real(RWP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RWP), intent(in)  :: V(:,:)         !< 1D eigenvalues
    real(RWP), intent(in)  :: W(:,:)         !< 1D weights
    real(RWP), intent(in)  :: g(:,:)         !< metrics
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs
    real(RWP), intent(in)  :: lambda         !< λ
    real(RWP), intent(in)  :: nu(:)          !< ν + νˢ
    real(RWP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RWP), intent(out) :: u(:,:,:,:)     !< solution

    integer :: nc, nd, np

    np = size(S,1)
    nc = size(S,3)
    nd = size(f,4)

    select case(np)
    case( 2)
      call TPO_Schwarz_I_Hand__2  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 3)
      call TPO_Schwarz_I_Hand__3  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 4)
      call TPO_Schwarz_I_Hand__4  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 5)
      call TPO_Schwarz_I_Hand__5  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 6)
      call TPO_Schwarz_I_Hand__6  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 7)
      call TPO_Schwarz_I_Hand__7  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
#ifndef __LIBXSMM__
    case( 8)
      call TPO_Schwarz_I_Hand__8  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case( 9)
      call TPO_Schwarz_I_Hand__9  (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(10)
      call TPO_Schwarz_I_Hand__10 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(11)
      call TPO_Schwarz_I_Hand__11 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(12)
      call TPO_Schwarz_I_Hand__12 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(13)
      call TPO_Schwarz_I_Hand__13 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(14)
      call TPO_Schwarz_I_Hand__14 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(15)
      call TPO_Schwarz_I_Hand__15 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(16)
      call TPO_Schwarz_I_Hand__16 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(17)
      call TPO_Schwarz_I_Hand__17 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(18)
      call TPO_Schwarz_I_Hand__18 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(19)
      call TPO_Schwarz_I_Hand__19 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(20)
      call TPO_Schwarz_I_Hand__20 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(21)
      call TPO_Schwarz_I_Hand__21 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(22)
      call TPO_Schwarz_I_Hand__22 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(23)
      call TPO_Schwarz_I_Hand__23 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(24)
      call TPO_Schwarz_I_Hand__24 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(25)
      call TPO_Schwarz_I_Hand__25 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(26)
      call TPO_Schwarz_I_Hand__26 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(27)
      call TPO_Schwarz_I_Hand__27 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(28)
      call TPO_Schwarz_I_Hand__28 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(29)
      call TPO_Schwarz_I_Hand__29 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(30)
      call TPO_Schwarz_I_Hand__30 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(31)
      call TPO_Schwarz_I_Hand__31 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(32)
      call TPO_Schwarz_I_Hand__32 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
    case(33)
      call TPO_Schwarz_I_Hand__33 (nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
#endif
    case default
      call TPO_Schwarz_I_Gen(S, V, W, g, cfg, lambda, nu, f, u)
    end select

  end subroutine TPO_Schwarz_I_Hand_RWP

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
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_

#ifndef __LIBXSMM__

#define _NP_  8
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__schwarz__3d_i__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxQt.F"
#include "IxQtxI.F"
#include "QtxIxI.F"

!===============================================================================

end module TPO__Schwarz__3D_I__Hand_RSP
