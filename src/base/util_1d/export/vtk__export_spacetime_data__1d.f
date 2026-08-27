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

!> summary:  Export 1D mesh space-time data into VTK XML file
!> author:   Erik Pfister
!> date:     2023/09/25
!===============================================================================

module VTK__Export_Spacetime_Data__1D
  use Kind_Parameters, only: RNP
  use Execution_Control
  use C_Binding
  use VTK_Binding
  use Structured_Mesh_Indexing__3D
  implicit none
  private

  public :: VTK_PackSpacetimeData_1D
  public :: VTK_ExportSpacetimeData_1D

contains

  !-----------------------------------------------------------------------------
  !> Packs variables of the same spacetime grid in one scalar array.
  !>
  !> This routine takes variables 'u' with corresponding names 'uname' and
  !> extends the arrays 's' and 'sname'. This routine can be called multiple
  !> times before the routine VTK_ExportSpacetimeData_1D is called.

  subroutine VTK_PackSpacetimeData_1D(u, s, uname, sname)

    real(RNP), intent(in) :: u(:,:,:,:,:)
      !< variables to be appended to the array `s` (0:po,ne,nc,0:pt,nt)
    real(RNP), allocatable, intent(inout) :: s(:,:,:,:,:)
      !< array of scalars to be extended (0:po,ne,0:pt,nt,ns->ns_new)
    character(len=*), intent(in) :: uname(:)
      !< names of the variables appended (nc)
    character(len=20), allocatable, intent(inout) :: sname(:)
      !< names of scalars contained in `s`(ns->ns_new)

    ! auxiliary variables ......................................................

    real(RNP)        , allocatable :: tmp_s(:,:,:,:,:)
    character(len=20), allocatable :: tmp_sname(:)

    integer :: nc     ! number of components appended
    integer :: ns     ! old number of scalars
    integer :: ns_new ! new number of scalars
    integer :: i

    ! determine number of components in u and current number of scalars in s
    nc = ubound(u,3)
    if (allocated(s)) then
      ns = size(s, 5)
    else
      ns = 0
    endif

    ! check if the length of uname matches the number of components
    if (size(uname) /= nc) then
      call Error( 'VTK_PackSpacetimeData_1D'                                       &
                , 'Length of uname does not match number of components in u' &
                , 'VTK__Export_Spacetime_Data__1D'                            )
    endif

    ! update ns
    ns_new = ns + nc

    ! check if s is not already allocated, allocate it with the correct dimensions
    if (.not. allocated(s)) then
      allocate(s(ubound(u,1),ubound(u,2),ubound(u,4),ubound(u,5),ns_new))
      allocate(sname(ns_new))
    else
      ! reallocate s to accommodate the new variables
      allocate(tmp_s(ubound(s,1),ubound(s,2),ubound(s,3),ubound(s,4),ns_new))
      allocate(tmp_sname(ns_new))
      tmp_s(:,:,:,:,1:ns) = s
      tmp_sname(1:ns)     = sname
      call move_alloc(tmp_s,     s    )
      call move_alloc(tmp_sname, sname)
    endif

    ! append new scalar variables for s
    do i = 1, nc
      s(:,:,:,:,ns+i) = u(:,:,i,:,:)
      sname(ns+i)     = uname(i)
    enddo

  end subroutine VTK_PackSpacetimeData_1D

  !-----------------------------------------------------------------------------
  !> Export elements and variables of mesh partition into VTK XML file.
  !>
  !> It is assumed that the mesh points are the element Lobatto nodes.
  !>
  !> This routine creates .vtu files of variables on a space-time grid
  !> for visualizations on ParaView. One axis represtens a spatial dimension
  !> and the other axis represents the evolution of that variable in time.

  subroutine VTK_ExportSpacetimeData_1D(x, t, s, sname, file)
    real(RNP), intent(in) :: x(0:,:)         !< mesh points in space (0:po,ne)
    real(RNP), intent(in) :: t(0:,:)         !< time points (0:pt,1:nt)
    real(RNP), intent(in) :: s(0:,:,0:,:,:)  !< scalars (0:po,ne,0:pt,nt,ns)
    character(len=*), intent(in) :: sname(:) !< scalar names (ns)
    character(len=*), intent(in) :: file     !< VTK output file

    ! VTK data .................................................................

    integer(C_INT) :: cell_type = VTK_QUAD

    integer(C_INT), allocatable :: cells(:,:)
    real(C_DOUBLE), allocatable :: points(:,:)
    real(C_DOUBLE), allocatable :: scalars(:,:)

    ! auxiliary variables ......................................................

    integer :: po, pt
    integer :: ne, np, ns, nt

    ! initialization ...........................................................

    po = ubound(x,1)
    ne = ubound(x,2)
    pt = ubound(t,1)
    nt = ubound(t,2)
    ns = size(s,5)

    ! points ...................................................................

    ! for a 2D space-time grid, calculate the number of spatial points and
    ! multiply by the number of time steps
    np = ne * (po + 1) * nt * (pt + 1)

    allocate(points(3,np))
    call BuildLinearSpacetimeCoords(np, po, ne, pt, nt, x, t, points)

    ! cells ....................................................................

    cell_type = VTK_QUAD
    call BuildLinearSpacetimeCells(po, pt, ne, nt, cells)

    ! scalars ..................................................................

    allocate(scalars(np,ns))
    call BuildLinearScalarData(np, ns, s, scalars)

    ! write ...................................................................

    call VTK_WriteXML_Unstructured( points, cells, cell_type &
                                  , ps       = scalars       &
                                  , ps_names = sname         &
                                  , file     = file          )

  end subroutine VTK_ExportSpacetimeData_1D

  !-----------------------------------------------------------------------------
  !> Maps spacetime element points to linear grid cells

  subroutine BuildLinearSpacetimeCoords(np, po, ne, pt, nt, xc, tc, points)
    integer,        intent(in)  :: np           !< number of VTK grid points
    integer,        intent(in)  :: po           !< polynomial degree in space
    integer,        intent(in)  :: ne           !< number of elements in space
    integer,        intent(in)  :: pt           !< polynomial degree in time
    integer,        intent(in)  :: nt           !< number of elements in time
    real(RNP),      intent(in)  :: xc(0:po,ne)  !< space mesh element points
    real(RNP),      intent(in)  :: tc(0:pt,nt)  !< time mesh element points
    real(C_DOUBLE), intent(out) :: points(3,np) !< VTK grid points

    integer :: i, j, k, l, idx

    idx = 1
    do l = 1, nt
      do k = 0, pt
        do j = 1, ne
          do i = 0, po
            points(1, idx) = xc(i, j)
            points(2, idx) = tc(k, l)
            points(3, idx) = 0
            idx = idx + 1
          end do
        end do
      end do
    end do

  end subroutine BuildLinearSpacetimeCoords

  !-----------------------------------------------------------------------------
  !> Connectivity of bilinear quadrilateral cells for spacetime grid.

  subroutine BuildLinearSpacetimeCells(po, pt, ne, nt, cells)
    integer, intent(in) :: po !< order of spatial elements
    integer, intent(in) :: pt !< order of temporal elements
    integer, intent(in) :: ne !< number of spatial elements
    integer, intent(in) :: nt !< number of temporal elements
    integer(C_INT), allocatable, intent(out) :: cells(:,:) !< grid cells

    integer :: c, o, i, j, l, m
    integer :: i0, j0, npe, npt

    allocate(cells(4, ne * nt * po * pt))

    npe = ne*(po+1) - 1
    npt = nt*(pt+1) - 1

    c =  1 ! cells counter
    o = -1 ! cells offset

    do m = 1, nt
      j0 = (m-1)*(pt+1)
      do j = 1, pt
      do l = 1, ne
        i0 = (l-1)*(po+1)
        do i = 1, po
          cells(1,c) = o + LexicalVertexIndex(i-1 + i0, j-1 + j0, 0, npe, npt, 0)
          cells(2,c) = o + LexicalVertexIndex(i   + i0, j-1 + j0, 0, npe, npt, 0)
          cells(3,c) = o + LexicalVertexIndex(i   + i0, j   + j0, 0, npe, npt, 0)
          cells(4,c) = o + LexicalVertexIndex(i-1 + i0, j   + j0, 0, npe, npt, 0)
          c = c + 1
        end do
      end do
      end do
    end do

  end subroutine BuildLinearSpacetimeCells

  !-----------------------------------------------------------------------------
  !> Maps scalar element variables to linear grid cells

  subroutine BuildLinearScalarData(np, ns, sc, scalars)
      real(RNP),      intent(in)  :: sc(:,:,:,:,:)
      real(C_DOUBLE), intent(out) :: scalars(:,:)
      integer,        intent(in)  :: np, ns

      scalars = reshape(sc, shape=[np, ns])

  end subroutine BuildLinearScalarData

  !=============================================================================

end module VTK__Export_Spacetime_Data__1D
