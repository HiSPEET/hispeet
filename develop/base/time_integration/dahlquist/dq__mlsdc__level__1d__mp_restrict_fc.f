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

!> summary:  Fine-to-coarse restriction of residual-like variables
!> author:   Erik Pfister, Jörg Stiller
!> date:     2024/03/22
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_Restrict_FC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse restriction of residual-like variables

  module subroutine Restrict_FC(this, r_f, r_c)
    class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
    complex(RNP), intent(in)    :: r_f(0:,:) !< fine residual variable
    complex(RNP), intent(inout) :: r_c(0:,:) !< coarse residual variable

    integer :: pt_c, pt_f ! polynomial degree in time
    integer :: mt_c, mt_f ! number of subintervals per time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: o_t        ! offset of time collocation points
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: l, m, n, n1, n2

    associate(iop_t => this % iop_cf_t)

      ! initialization .........................................................

      ! coarse dimensions
      mt_c = ubound(r_c, 1)
      nt_c = ubound(r_c, 2)

      ! fine dimensions
      mt_f = ubound(r_f, 1)
      nt_f = ubound(r_f, 2)

      ! polynomial degree in time and offset collocation points
      select case(this % sdc % nodes)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        o_t  = 1
      case('E','L')
        ! equidistant of Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        o_t  = 0
      end select

      ! consistency of temporal dimensions
      is_consistent = pt_c == iop_t % po_c .and. &
                      pt_f == iop_t % po_f

      ! consistency with spatial interpolation mode
      select case(iop_t % mode)
      case(1)
        is_consistent = is_consistent .and. nt_f == nt_c
      case(2)
        is_consistent = is_consistent .and. nt_f == nt_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Restrict_FC'              &
                  , 'failed consistency check' &
                  , 'DQ__MLSDC__Level__1D'     )
      end if

      ! temporal interpolation .................................................

      select case(iop_t % mode)

      case(0)
        
        r_c = r_f

      case(1)
        
        ! nt_f = nt_c
        if (o_t == 1) then
          r_c(0,:) = (0.0, 0.0)
        end if
        do n = 1  , nt_c
        do m = o_t, mt_c
          r_c(m,n) = (0.0, 0.0)
          !r_c(:,n)%re = matmul(transpose(iop_t%A(:,:,1)),r_f(:,n)%re)
          !r_c(:,n)%im = matmul(transpose(iop_t%A(:,:,1)),r_f(:,n)%im)
          do l = o_t, mt_f
            r_c(m,n)%re = r_c(m,n)%re + iop_t % A(l-o_t,m-o_t,1) * r_f(l,n)%re
            r_c(m,n)%im = r_c(m,n)%im + iop_t % A(l-o_t,m-o_t,1) * r_f(l,n)%im
          end do
        end do
        end do

      case(2)

      ! nt_f = 2 * nt_c

      ! maybe TBD
      print*, "Not supported yet."
      stop

      end select

    end associate

  end subroutine Restrict_FC

  !=============================================================================

end submodule MP_Restrict_FC
