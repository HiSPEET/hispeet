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

!> summary:  Fine-to-coarse space-time solution projection
!> author:   Erik Pfister, Joerg Stiller
!> date:     2024/03/22
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_Project_FC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse space-time projection of solution-like variables

  module subroutine Project_FC(this, u_f, u_c)
    class(DQ_MLSDC_Level_1D), intent(in) :: this !< coarse level
    complex(RNP), intent(in)    :: u_f(0:,:) !< fine solution variable
    complex(RNP), intent(inout) :: u_c(0:,:) !< coarse solution variable

    integer :: pt_c, pt_f ! polynomial degree in time
    integer :: mt_c, mt_f ! number of subintervals per time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: o_t        ! offset of time collocation points
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: l, m, n

    associate(pop_t => this % pop_fc_t)

      ! initialization .........................................................

      ! fine dimensions
      mt_f = ubound(u_f, 1)
      nt_f = ubound(u_f, 2)

      ! coarse dimensions
      mt_c = ubound(u_c, 1)
      nt_c = ubound(u_c, 2)

      ! polynomial degree in time and offset collocation points
      select case(this % sdc % nodes)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        o_t   = 1
      case('E','L')
        ! equidistant or Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        o_t   = 0
      end select

      ! consistency of temporal dimensions
      is_consistent = pt_f == pop_t % po_f .and. &
                      pt_c == pop_t % po_c

      ! consistency with spatial projection mode
      select case(pop_t % mode)
      case(1)
        is_consistent = is_consistent .and. nt_f == nt_c
      case(2)
        is_consistent = is_consistent .and. nt_f == nt_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Project_FC'               &
                  , 'failed consistency check' &
                  , 'DQ__MLSDC__Level__1D'     )
      end if

      ! temporal projection ....................................................

      select case(pop_t % mode)

      case(0)
        
        u_c = u_f

      case(1)
        ! nt_f = nt_c
        do n = 1, nt_c

          ! inject boundary values (m = 0, mt)
          u_c(0    ,n) = u_f(0    ,n)
          u_c(mt_c ,n) = u_f(mt_f ,n)

          do m = 1, mt_c-1
          u_c(m,n) = (0.0, 0.0)
          do l = o_t, mt_f
            u_c(m,n)%re = u_c(m,n)%re + pop_t % A(m-o_t,l-o_t,1) * u_f(l,n)%re
            u_c(m,n)%im = u_c(m,n)%im + pop_t % A(m-o_t,l-o_t,1) * u_f(l,n)%im
          end do
          end do
        end do

      case(2)

        ! nt_f = 2 * nt_c

        ! maybe TBD
        print*, "Not supported."
        stop

      end select

    end associate

  end subroutine Project_FC

  !=============================================================================

end submodule MP_Project_FC
