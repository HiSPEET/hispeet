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

!> summary:  Coarse-to-fine space-time solution interpolation
!> author:   Erik Pfister, Joerg Stiller
!> date:     2023/08/25
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_Interpolate_CF
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Coarse-to-fine space-time interpolation of solution-like variables

  module subroutine Interpolate_CF(this, u_c, u_f, complete)
    class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
    complex(RNP), intent(in)    :: u_c(0:,:) !< coarse solution variable
    complex(RNP), intent(inout) :: u_f(0:,:) !< fine solution variable
    logical,      intent(in)    :: complete  !< T/F for all/refined elements

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
      mt_c = ubound(u_c, 1)
      nt_c = ubound(u_c, 2)

      ! fine dimensions
      mt_f = ubound(u_f, 1)
      nt_f = ubound(u_f, 2)

      ! polynomial degree in time and offset collocation points
      select case(this % sdc % nodes)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        o_t  = 1
      case('E','L')
        ! equidistant or Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        o_t  = 0
      end select

      ! consistency of temporal dimensions
      is_consistent = pt_c == iop_t % po_c .and. &
                      pt_f == iop_t % po_f

      ! consistency with temporal interpolation mode
      select case(iop_t % mode)
      case(1)
        is_consistent = is_consistent .and. nt_f == nt_c
      case(2)
        is_consistent = is_consistent .and. nt_f == nt_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Interpolate_CF'           &
                  , 'failed consistency check' &
                  , 'DQ__MLSDC__Level__1D'     )
      end if

      ! temporal interpolation .................................................

      select case(iop_t % mode)

      case(0)

        u_f = u_c

      case(1)

        ! nt_f = nt_c
        do n = 1, nt_c

          ! inject boundary values (m = 0, mt)
          ! u(0,nt) = u(mt,nt-1)
          ! The first sub-time point cannot be changed; Dahlquist has one step.
          ! Dahlquist IC
          u_f(0    ,n) = u_c(0    ,n)
          u_f(mt_f ,n) = u_c(mt_c ,n)

          ! interpolate interior values
          do m = 1, mt_f-1
          u_f(m,n) = (0.0, 0.0)
          do l = 0, pt_c
            !! offset o_t = 1 required with right-sided points !!
            u_f(m,n)%re = u_f(m,n)%re + iop_t % A(m-o_t,l,1) * u_c(l+o_t,n)%re
            u_f(m,n)%im = u_f(m,n)%im + iop_t % A(m-o_t,l,1) * u_c(l+o_t,n)%im
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

  end subroutine Interpolate_CF

  !=============================================================================

end submodule MP_Interpolate_CF
