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

!> summary:  3D multilevel spacetime mesh variable
!> author:   Joerg Stiller
!> date:     2024/11/24
!===============================================================================

module ML__Spacetime_Variable__3D
  use Mesh_Variable__3D
  use Spacetime_Variable__3D
  use ML__Spacetime_Operators__3D
  implicit none
  private

  public :: ML_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spacetime mesh variable
  !>
  !> The values of a multilevel spacetime variable `u` are accessed via
  !>
  !>       u % level(l) % var(m,n) % val(i,j,k,e,c)
  !>
  !>  where
  !>
  !>    - `1 ≤   l   ≤ l_top`:  level
  !>    - `0 ≤   m   ≤ pt(l)`:  temporal collocation point ID
  !>    - `1 ≤   n   ≤ nt(l)`:  temporal element ID in given time slice
  !>    - `0 ≤ i,j,k ≤ po(l)`:  spatial collocation point triple index
  !>    - `1 ≤   e   ≤ ne(l)`:  spatial element ID
  !>    - `1 ≤   c   ≤ nc   `:  component ID

  type ML_SpacetimeVariable_3D
    integer :: nc !< number of components
    type(SpacetimeVariable_3D), allocatable :: level(:) !< variable per level
    character(len=:),           allocatable :: name(:)  !< component names
  contains
    procedure :: Init => Init_ML_SpacetimeVariable_3D
    procedure :: GetSlice
  end type ML_SpacetimeVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel spacetime variable

  subroutine Init_ML_SpacetimeVariable_3D(this, ml_op, nc, name)
    class(ML_SpacetimeVariable_3D),  intent(inout) :: this
    class(ML_SpacetimeOperators_3D), intent(in)    :: ml_op
    integer,                         intent(in)    :: nc
    character(len=*),      optional, intent(in)    :: name(nc)

    integer :: c, l

    allocate(this % level( size(ml_op%sem) ))

    this % nc = nc

    do l = 1, size(this%level)
      call this % level(l) % Init( mesh = ml_op % sem(l) % mesh        &
                                 , po_x = ml_op % sem(l) % std_op % po &
                                 , po_t = ml_op % po_t(l)              &
                                 , ne_t = ml_op % ne_t(l)              &
                                 , nc   = nc                           )
    end do

    if (present(name)) then
      this % name = name
    else
      l = 2 + int(log10(dble(nc)))
      allocate(character(len=l) :: this % name(nc))
      do c = 1, nc
        write(this%name(c), '(A,I0)') 'v', c
      end do
    end if

  end subroutine Init_ML_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> Create a new multilevel spacetime variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % var(m,n) % val` if `copy`
  !> is false or absent. Otherwise they are stored in fresh memory.

  subroutine GetSlice(this, slice, first, last, copy)
    class(ML_SpacetimeVariable_3D), intent(in)    :: this
    class(ML_SpacetimeVariable_3D), intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    integer :: l

    if (allocated(slice % level)) then
      deallocate(slice % level)
    end if

    allocate(slice % level(size(this%level)))

    slice % nc   = last - first + 1
    slice % name = this % name(first:last)

    do l = 1, size(this%level)
      call this % level(l) % GetSlice(slice%level(l), first, last, copy)
    end do

  end subroutine GetSlice

  !=============================================================================

end module ML__Spacetime_Variable__3D
