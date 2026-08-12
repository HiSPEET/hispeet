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

!> summary:  3D spacetime interpolation from coarse to fine level
!> author:   Joerg Stiller
!> date:     2024/11/29
!===============================================================================

module CF__Spacetime_Interpolation__3D
  use Constants, only: ZERO
  use Array_Assignments
  use Mesh_Variable__3D
  use Spacetime_Variable__3D
  use Parent_To_Child_Interpolation__3D
  use ML__Spacetime_Operators__3D

  implicit none
  private

  public :: CF_SpacetimeInterpolation_3D

contains

  !-----------------------------------------------------------------------------
  !> Coarse-to-fine space-time interpolation of solution-like variables

  subroutine CF_SpacetimeInterpolation_3D(ml_op, l_c, v_c, v_f, skip_frozen)
    class(ML_SpacetimeOperators_3D), intent(in)    :: ml_op
    integer,                         intent(in)    :: l_c !< coarse level ID
    class(SpacetimeVariable_3D),     intent(in)    :: v_c !< coarse variable
    class(SpacetimeVariable_3D),     intent(inout) :: v_f !< fine variable
    logical, optional, intent(in) :: skip_frozen !< exclude frozen elements [F]

    type(SpacetimeVariable_3D), allocatable, save :: v_i
    logical :: skip_frozen_
    integer :: l_f
    integer :: px_c, pt_c, nt_c
    integer :: px_f, nx_f, pt_f, nt_f
    integer :: c, e, k, m, n, n1, n2, nc

    ! initialization ...........................................................

    if (present(skip_frozen)) then
      skip_frozen_ = skip_frozen
    else
      skip_frozen_ = .false.
    end if

    l_f = l_c + 1

    ! coarse mesh dimensions
    px_c = ml_op % sem(l_c) % std_op % po
    pt_c = v_c % po_t
    nt_c = v_c % ne_t

    ! fine mesh dimensions
    px_f = ml_op % sem(l_f) % std_op % po
    if (skip_frozen_) then
      nx_f = ml_op % sem(l_f) % mesh % n_elem_active
    else
      nx_f = ml_op % sem(l_f) % mesh % n_elem
    end if
    pt_f = v_f % po_t
    nt_f = v_f % ne_t

    ! number of components
    nc = v_c % nc

    !$omp master
    allocate(v_i)
    v_i % po_t = pt_c
    v_i % ne_t = nt_c
    allocate(v_i % var(0:pt_c, 1:nt_c))
    do n = 1, nt_c
    do m = 0, pt_c
      call v_i % var(m,n) % Init(ml_op % sem(l_f) % mesh, px_f, nc)
    end do
    end do
    !$omp end master

    ! spatial interpolation ....................................................

    do n = 1, nt_c
    do m = 0, pt_c
      call ParentToChildInterpolation_3D( ml_op % sem(l_c) % mesh & ! parent
                                        , ml_op % sem(l_f) % mesh & ! child
                                        , ml_op % iop_cf_x(l_c)   & ! iop
                                        , v_c % var(m,n) % val    & ! v_p
                                        , v_i % var(m,n) % val    & ! v_c
                                        , skip_frozen_            )
    end do
    end do

    ! temporal interpolation ...................................................

    associate(A => ml_op % iop_cf_t(l_c) % A)

      select case(size(A,3))

      case(0)

        ! identity: pt_f = pt_c, nt_f = nt_c
        do n = 1, nt_c
        do m = 0, pt_c
          do c = 1, nc
            call SetArray( v_f % var(m,n) % val(:,:,:,:nx_f,c) &
                         , v_i % var(m,n) % val(:,:,:,:nx_f,c) )
          end do
        end do
        end do

      case(1)

        ! p-refinement: nt_f = nt_c
        do n = 1, nt_c
        do m = 0, pt_f
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, nx_f
            v_f % var(m,n) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_c
              v_f % var(m,n) % val(:,:,:,e,c)                   &
                  = v_f % var(m,n) % val(:,:,:,e,c)             &
                  + v_i % var(k,n) % val(:,:,:,e,c) * A(m,k,1)
            end do
          end do
          end do
        end do
        end do

      case(2)

        ! hp-refinement: nt_f = 2 * nt_c
        do n = 1, nt_c
        do m = 0, pt_f
          n1 = 2 * n - 1
          n2 = 2 * n
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, nx_f
            v_f % var(m,n1) % val(:,:,:,e,c) = ZERO
            v_f % var(m,n2) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_c
              v_f % var(m,n1) % val(:,:,:,e,c)                   &
                  = v_f % var(m,n1) % val(:,:,:,e,c)             &
                  + v_i % var(k,n ) % val(:,:,:,e,c) * A(m,k,1)
              v_f % var(m,n2) % val(:,:,:,e,c)                   &
                  = v_f % var(m,n2) % val(:,:,:,e,c)             &
                  + v_i % var(k,n ) % val(:,:,:,e,c) * A(m,k,2)
            end do
          end do
          end do
        end do
        end do

      end select
    end associate

    ! finalization .............................................................

    !$omp master
    deallocate(v_i)
    !$omp end master

  end subroutine CF_SpacetimeInterpolation_3D

  !=============================================================================

end module CF__Spacetime_Interpolation__3D
