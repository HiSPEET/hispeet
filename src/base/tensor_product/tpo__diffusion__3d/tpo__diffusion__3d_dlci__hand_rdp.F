!> summary:   element diffusion operator based on handcrafted kernels (DLCI)
!> author:    Jerome Michel, Joerg Stiller
!> date:      2021/10/12
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_DLCI__Hand_RDP
  use Kind_Parameters, only: RWP => RDP
  use TPO__Diffusion__3D_DLCI__Gen
  implicit none
  private

  public :: TPO_Diffusion_DLCI_Hand

  interface TPO_Diffusion_DLCI_Hand
    module procedure TPO_Diffusion_DLCI_Hand_RWP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> element diffusion operator based on handcrafted kernels (DLCI)

  subroutine TPO_Diffusion_DLCI_Hand_RWP(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                         Ji_n, ub, qb)

    real(RWP), intent(in)  :: Ms(:)        !< 1D standard mass matrix       (np)
    real(RWP), intent(in)  :: Ds(:,:)      !< 1D standard diff matrix    (np,np)
    real(RWP), intent(in)  :: Jd(:,:,:,:)  !< Jacobian determinant (np,np,np,ne)
    real(RWP), intent(in)  :: G(:,:,:,:,:) !< Laplacian metrics
    real(RWP), intent(in)  :: lambda       !< Helmholtz parameter
    real(RWP), intent(in)  :: nu           !< diffusivity
    real(RWP), intent(in)  :: u(:,:,:,:)   !< operand              (np,np,np,ne)
    real(RWP), intent(out) :: v(:,:,:,:)   !< result               (np,np,np,ne)

    real(RWP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
    real(RWP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RWP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    integer :: np, ne

    np = size(Ms)
    ne = size(u,4)

    select case(np)
    case( 2)
      call TPO_Diffusion_DLCI_Hand__2  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 3)
      call TPO_Diffusion_DLCI_Hand__3  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 4)
      call TPO_Diffusion_DLCI_Hand__4  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 5)
      call TPO_Diffusion_DLCI_Hand__5  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 6)
      call TPO_Diffusion_DLCI_Hand__6  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 7)
      call TPO_Diffusion_DLCI_Hand__7  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
#ifndef __LIBXSMM__
    case( 8)
      call TPO_Diffusion_DLCI_Hand__8  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case( 9)
      call TPO_Diffusion_DLCI_Hand__9  (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(10)
      call TPO_Diffusion_DLCI_Hand__10 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(11)
      call TPO_Diffusion_DLCI_Hand__11 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(12)
      call TPO_Diffusion_DLCI_Hand__12 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(13)
      call TPO_Diffusion_DLCI_Hand__13 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(14)
      call TPO_Diffusion_DLCI_Hand__14 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(15)
      call TPO_Diffusion_DLCI_Hand__15 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(16)
      call TPO_Diffusion_DLCI_Hand__16 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(17)
      call TPO_Diffusion_DLCI_Hand__17 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(18)
      call TPO_Diffusion_DLCI_Hand__18 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(19)
      call TPO_Diffusion_DLCI_Hand__19 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(20)
      call TPO_Diffusion_DLCI_Hand__20 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(21)
      call TPO_Diffusion_DLCI_Hand__21 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(22)
      call TPO_Diffusion_DLCI_Hand__22 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(23)
      call TPO_Diffusion_DLCI_Hand__23 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(24)
      call TPO_Diffusion_DLCI_Hand__24 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(25)
      call TPO_Diffusion_DLCI_Hand__25 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(26)
      call TPO_Diffusion_DLCI_Hand__26 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(27)
      call TPO_Diffusion_DLCI_Hand__27 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(28)
      call TPO_Diffusion_DLCI_Hand__28 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(29)
      call TPO_Diffusion_DLCI_Hand__29 (ne, Ms, Ds, Jd, G, lambda, nu, u, v), &
                                        Ji_n, ub, qb)
    case(30)
      call TPO_Diffusion_DLCI_Hand__30 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(31)
      call TPO_Diffusion_DLCI_Hand__31 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(32)
      call TPO_Diffusion_DLCI_Hand__32 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
    case(33)
      call TPO_Diffusion_DLCI_Hand__33 (ne, Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)
#endif
    case default
      call TPO_Diffusion_DLCI_Gen(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                  Ji_n, ub, qb)
    end select

  end subroutine TPO_Diffusion_DLCI_Hand_RWP

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
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  3
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  4
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  5
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  6
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  7
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_

#ifndef __LIBXSMM__

#define _NP_  8
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_  9
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 10
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 11
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 12
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 13
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 14
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 15
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 16
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 17
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 18
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 19
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 20
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 21
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 22
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 23
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 24
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 25
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 26
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 27
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 28
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 29
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 30
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 31
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 32
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_
#define _NP_ 33
#include "tpo__diffusion__3d_dlci__hand_kernel.f"
#undef _NP_

#endif

!===============================================================================
! primitives for direction 1,2,3

#include "IxIxQt.F"
#include "IxQtxI.F"
#include "QtxIxI.F"

!===============================================================================

end module TPO__Diffusion__3D_DLCI__Hand_RDP
