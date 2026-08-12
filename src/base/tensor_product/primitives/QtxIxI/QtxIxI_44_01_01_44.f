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

!-----------------------------------------------------------------------------
!> Computes  v = α Qᵀ⊗I⊗I u + β v
!>
!>   * k:  block 4,  unroll 4
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 4,  unroll 4
!>   * explicit remainder handling (r)

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(QtxIxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_*_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_*_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p   ! loop counters
  integer   :: kb, pb       ! block counters

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do ij = 1, _NQ_**2
    v(ij,k) = beta * v(ij,k)
  end do
  end do

#if _NQ_T4_ > 0

  !$acc loop independent vector
  do kb = 1, NQ_T4, 4
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop collapse(2) independent vector
      do ij  = 1, _NQ_**2

        tmp0 =        Q(pb  ,kb  ) * u(ij,pb  )
        tmp1 =        Q(pb  ,kb+1) * u(ij,pb  )
        tmp2 =        Q(pb  ,kb+2) * u(ij,pb  )
        tmp3 =        Q(pb  ,kb+3) * u(ij,pb  )

        tmp0 = tmp0 + Q(pb+1,kb  ) * u(ij,pb+1)
        tmp1 = tmp1 + Q(pb+1,kb+1) * u(ij,pb+1)
        tmp2 = tmp2 + Q(pb+1,kb+2) * u(ij,pb+1)
        tmp3 = tmp3 + Q(pb+1,kb+3) * u(ij,pb+1)

        tmp0 = tmp0 + Q(pb+2,kb  ) * u(ij,pb+2)
        tmp1 = tmp1 + Q(pb+2,kb+1) * u(ij,pb+2)
        tmp2 = tmp2 + Q(pb+2,kb+2) * u(ij,pb+2)
        tmp3 = tmp3 + Q(pb+2,kb+3) * u(ij,pb+2)

        tmp0 = tmp0 + Q(pb+3,kb  ) * u(ij,pb+3)
        tmp1 = tmp1 + Q(pb+3,kb+1) * u(ij,pb+3)
        tmp2 = tmp2 + Q(pb+3,kb+2) * u(ij,pb+3)
        tmp3 = tmp3 + Q(pb+3,kb+3) * u(ij,pb+3)

        v(ij,kb  ) = alpha * tmp0 + v(ij,kb  )
        v(ij,kb+1) = alpha * tmp1 + v(ij,kb+1)
        v(ij,kb+2) = alpha * tmp2 + v(ij,kb+2)
        v(ij,kb+3) = alpha * tmp3 + v(ij,kb+3)

      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! 1:nq/4*4, nq
  ! -----------------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do kb  = 1, NQ_T4, 4
    do ij   = 1, _NQ_**2

      v(ij,kb  ) = alpha * Q(pb,kb  ) * u(ij,pb) + v(ij,kb  )
      v(ij,kb+1) = alpha * Q(pb,kb+1) * u(ij,pb) + v(ij,kb+1)
      v(ij,kb+2) = alpha * Q(pb,kb+2) * u(ij,pb) + v(ij,kb+2)
      v(ij,kb+3) = alpha * Q(pb,kb+3) * u(ij,pb) + v(ij,kb+3)

    end do
  end do

  ! nq, 1:nq
  ! ---------------------------------------------------------------

  kb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do ij = 1, _NQ_**2
    do p = 1, _NQ_

        v(ij,kb) = alpha * Q(p,kb) * u(ij,p) + v(ij,k)

    end do
  end do


#elif _NQ_ == _NQ_T4_ + 2

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do kb = 1, NQ_T4, 4
    do ij  = 1, _NQ_**2

      tmp0 =        Q(pb  ,kb  ) * u(ij,pb  )
      tmp1 =        Q(pb  ,kb+1) * u(ij,pb  )
      tmp2 =        Q(pb  ,kb+2) * u(ij,pb  )
      tmp3 =        Q(pb  ,kb+3) * u(ij,pb  )

      tmp0 = tmp0 + Q(pb+1,kb  ) * u(ij,pb+1)
      tmp1 = tmp1 + Q(pb+1,kb+1) * u(ij,pb+1)
      tmp2 = tmp2 + Q(pb+1,kb+2) * u(ij,pb+1)
      tmp3 = tmp3 + Q(pb+1,kb+3) * u(ij,pb+1)

      v(ij,kb  ) = alpha * tmp0 + v(ij,kb  )
      v(ij,kb+1) = alpha * tmp1 + v(ij,kb+1)
      v(ij,kb+2) = alpha * tmp2 + v(ij,kb+2)
      v(ij,kb+3) = alpha * tmp3 + v(ij,kb+3)

    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  kb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = kb, kb+1
    do ij = 1, _NQ_**2
      do p = 1, _NQ_

        v(ij,k) = alpha * Q(p,k) * u(ij,p) + v(ij,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 3

 ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do kb = 1, NQ_T4, 4
    do ij  = 1, _NQ_**2

      tmp0 =        Q(pb  ,kb  ) * u(ij,pb  )
      tmp1 =        Q(pb  ,kb+1) * u(ij,pb  )
      tmp2 =        Q(pb  ,kb+2) * u(ij,pb  )
      tmp3 =        Q(pb  ,kb+3) * u(ij,pb  )

      tmp0 = tmp0 + Q(pb+1,kb  ) * u(ij,pb+1)
      tmp1 = tmp1 + Q(pb+1,kb+1) * u(ij,pb+1)
      tmp2 = tmp2 + Q(pb+1,kb+2) * u(ij,pb+1)
      tmp3 = tmp3 + Q(pb+1,kb+3) * u(ij,pb+1)

      tmp0 = tmp0 + Q(pb+2,kb  ) * u(ij,pb+2)
      tmp1 = tmp1 + Q(pb+2,kb+1) * u(ij,pb+2)
      tmp2 = tmp2 + Q(pb+2,kb+2) * u(ij,pb+2)
      tmp3 = tmp3 + Q(pb+2,kb+3) * u(ij,pb+2)

      v(ij,kb  ) = alpha * tmp0 + v(ij,kb  )
      v(ij,kb+1) = alpha * tmp1 + v(ij,kb+1)
      v(ij,kb+2) = alpha * tmp2 + v(ij,kb+2)
      v(ij,kb+3) = alpha * tmp3 + v(ij,kb+3)

    end do
  end do

  ! nq-2:nq, 1:nq ---------------

  kb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = kb, kb+2
    do ij = 1, _NQ_**2
      do p = 1, _NQ_

        v(ij,k) = alpha * Q(p,k) * u(ij,p) + v(ij,k)

      end do
    end do
  end do

#endif

end subroutine PROC(QtxIxI__,_NQ_)

#undef _NQ_T4
