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

!> summary:  Export 3d mesh data into VTK XML file
!> author:   Joerg Stiller
!> date:     2014/11/27, 2026/08/26
!===============================================================================

module VTK__Export_Mesh_Data__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Gauss_Jacobi,    only: LobattoPoints, LobattoPolynomial
  use VTK_Binding
  use TPO__AAA__3D
  use Structured_Mesh_Indexing__3D
  implicit none
  private

  public :: VTK_ExportMeshData_3D

contains

  !-----------------------------------------------------------------------------
  !> Export elements and variables of mesh partition into VTK XML file.
  !>
  !> It is assumed that the mesh points are the element Lobatto nodes.
  !>
  !> In case of partitioned (multi-piece) data, the partition ID, part, and the
  !> number of partitions, n_parts, must be given. The partitions are numbered
  !> from 0 to n_parts-1. Partition 0 writes the parallel unstructured mesh file
  !> (PVTU).

  subroutine VTK_ExportMeshData_3D( x, s, sname, v, vname, file &
                                    , part, n_parts, subdiv, mask )

    real(RNP), intent(in) :: x(0:,0:,0:,:,:)
    !< element nodes (0:po,0:po,0:po,ne,3)

    real(RNP), optional, intent(in) :: s(0:,0:,0:,:,:)
    !< scalars (0:po,0:po,0:po,ne,ns)
    character(len=*), optional, intent(in) :: sname(:)
    !< scalar names (ns)

    real(RNP), optional, intent(in) :: v(0:,0:,0:,:,:,:)
    !< vectors (0:po,0:po,0:po,ne,3,nv)
    character(len=*), optional, intent(in) :: vname(:)
    !< vector names (nv)

    ! output control
    character(len=*),  intent(in) :: file    !< VTK output file
    integer, optional, intent(in) :: part    !< partition (piece)
    integer, optional, intent(in) :: n_parts !< number of partitions (pieces)
    logical, optional, intent(in) :: subdiv  !< T: quadratic subdivision [auto]
    logical, optional, intent(in) :: mask(:) !< skip elements with mask(e) = F

    ! VTK data .................................................................

    integer(VTK_INT32) :: cell_type
    integer(VTK_INT32), allocatable :: cells(:,:)
    real(VTK_FLOAT64),  allocatable :: points(:,:)
    real(VTK_FLOAT64),  allocatable :: ps(:,:)
    real(VTK_FLOAT64),  allocatable :: pv(:,:,:)

    ! auxiliary variables ......................................................

    logical,   allocatable :: mask_(:) ! element mask
    real(RNP), allocatable :: iop(:,:) ! interpolation operator

    integer :: ep, ne, ni, np, ns, nv, nx, po
    integer :: interpolation_order

    ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(part)) then
      if (part < 0) return ! skip empty partition
    end if

    ! dimensions ...............................................................

    po = ubound(x,1)
    ne = ubound(x,4)
    ep = po + 1

    if (present(mask)) then
      allocate(mask_(ne), source = mask)
    else
      allocate(mask_(ne), source = .true.)
    end if

    ! number of exported elements
    nx = count(mask_)

    ! interpolation order ......................................................

    if (po < 1) then
      return
    else if (po == 1) then
      interpolation_order = 1
    else if (present(subdiv)) then
      if (subdiv) then
        interpolation_order = 2
      else
        interpolation_order = 1
      end if
    else
      interpolation_order = 2
    end if

    ! interpolation points per element and direction
    ni = interpolation_order * po + 1

    ! grid points ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    np = nx * ni**3
    allocate(points(3,np))

    select case(interpolation_order)
    case(1)
      call BuildLinearPointCoords(ep, ne, nx, mask_, x, points)
    case(2)
      call BuildInterpolationOperator(po, iop)
      call BuildQuadraticPointCoords(ep, ni, ne, nx, mask_, iop, x, points)
    end select

    ! mesh cells :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    select case(interpolation_order)
    case(1)
      cell_type = VTK_HEXAHEDRON
      call BuildLinearCells(po, nx, cells)
    case(2)
      cell_type = VTK_TRIQUADRATIC_HEXAHEDRON
      call BuildQuadraticCells(po, nx, cells)
    end select

    ! scalars ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(s) .and. present(sname)) then
      ns = size(sname)
      allocate(ps(np,ns))

      select case(interpolation_order)
      case(1)
        call BuildLinearScalarData(ep, ne, nx, ns, mask_, s, ps)
      case(2)
        call BuildQuadraticScalarData(ep, ni, ne, nx, ns, mask_, iop, s, ps)
      end select

    end if

    ! vectors ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(v) .and. present(vname)) then
      nv = size(vname)
      allocate(pv(3,np,nv))

      select case(interpolation_order)
      case(1)
        call BuildLinearVectorData(ep, ne, nx, nv, mask_, v, pv)
      case(2)
        call BuildQuadraticVectorData(ep, ni, ne, nx, nv, mask_, iop, v, pv)
      end select

    end if

    ! Export :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    call VTK_WriteXML_Unstructured( points, cells, cell_type  &
                                  , ps       = ps             &
                                  , ps_names = sname          &
                                  , pv       = pv             &
                                  , pv_names = vname          &
                                  , file     = file           &
                                  , piece    = part           &
                                  , n_pieces = n_parts        )

  end subroutine VTK_ExportMeshData_3D

  !-----------------------------------------------------------------------------
  !> Maps element points to linear grid cells

  subroutine BuildLinearPointCoords(ep, ne, nx, mask, xe, points)
    integer,   intent(in)  :: ep       !< element points per direction
    integer,   intent(in)  :: ne       !< num elements
    integer,   intent(in)  :: nx       !< num exported elements
    logical,   intent(in)  :: mask(ne) !< element mask
    real(RNP), intent(in)  :: xe(ep,ep,ep,ne,3) !< mesh points
    real(VTK_FLOAT64), intent(out) :: points(3,ep,ep,ep,nx) !< VTK mesh points

    integer :: e, i, j, k, l

    l = 0
    do e = 1, ne
      if (mask(e)) then
        l = l + 1
        do k = 1, ep
        do j = 1, ep
        do i = 1, ep
          points(1,i,j,k,l) = real(xe(i,j,k,e,1), VTK_FLOAT64)
          points(2,i,j,k,l) = real(xe(i,j,k,e,2), VTK_FLOAT64)
          points(3,i,j,k,l) = real(xe(i,j,k,e,3), VTK_FLOAT64)
        end do
        end do
        end do
      end if
    end do

  end subroutine BuildLinearPointCoords

  !-----------------------------------------------------------------------------
  !> Connectivity of trilinear hexahedral cells.

  subroutine BuildLinearCells(po, nx, cells)
    integer, intent(in) :: po  !< order of elements
    integer, intent(in) :: nx  !< number of exported elements
    integer(VTK_INT32), allocatable, intent(out) :: cells(:,:) !< VTK cells

    integer :: c, i, j, k, l, n, o

    allocate(cells(8, nx * po**3))

    n = (po + 1)**3  ! grid points per element
    o = -1           ! offset of point IDs
    c =  1           ! cell counter

    do l = 1, nx
      do k = 1, po
      do j = 1, po
      do i = 1, po
        cells(1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, po, po)
        cells(2,c) = o + LexicalVertexIndex(i  , j-1, k-1, po, po)
        cells(3,c) = o + LexicalVertexIndex(i  , j  , k-1, po, po)
        cells(4,c) = o + LexicalVertexIndex(i-1, j  , k-1, po, po)
        cells(5,c) = o + LexicalVertexIndex(i-1, j-1, k  , po, po)
        cells(6,c) = o + LexicalVertexIndex(i  , j-1, k  , po, po)
        cells(7,c) = o + LexicalVertexIndex(i  , j  , k  , po, po)
        cells(8,c) = o + LexicalVertexIndex(i-1, j  , k  , po, po)
        c = c + 1
      end do
      end do
      end do
      o = o + n
    end do

  end subroutine BuildLinearCells

  !-----------------------------------------------------------------------------
  !> Maps scalar element variables to linear grid cells

  subroutine BuildLinearScalarData(ep, ne, nx, ns, mask, se, ps)
    integer,   intent(in)  :: ep       !< element points per direction
    integer,   intent(in)  :: ne       !< num elements
    integer,   intent(in)  :: nx       !< num exported elements
    integer,   intent(in)  :: ns       !< num scalars
    logical,   intent(in)  :: mask(ne) !< element mask
    real(RNP), intent(in)  :: se(ep,ep,ep,ne,ns) !< SEM scalars
    real(VTK_FLOAT64), intent(out) :: ps(ep,ep,ep,nx,ns) !< VTK scalars

    integer :: e, i, j, k, l, n

    do n = 1, ns
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, ep
          do j = 1, ep
          do i = 1, ep
            ps(i,j,k,l,n) = real(se(i,j,k,e,n), VTK_FLOAT64)
          end do
          end do
          end do
        end if
      end do
    end do

  end subroutine BuildLinearScalarData

  !-----------------------------------------------------------------------------
  !> Maps vector element variables to linear grid cells

  subroutine BuildLinearVectorData(ep, ne, nx, nv, mask, ve, pv)
    integer,   intent(in)  :: ep       !< element points per direction
    integer,   intent(in)  :: ne       !< num elements
    integer,   intent(in)  :: nx       !< num exported elements
    integer,   intent(in)  :: nv       !< num vectors
    logical,   intent(in)  :: mask(ne) !< element mask
    real(RNP), intent(in)  :: ve(ep,ep,ep,ne,3,nv) !< SEM vectors
    real(VTK_FLOAT64), intent(out) :: pv(3,ep,ep,ep,nx,nv) !< VTK vectors

    integer :: e, i, j, k, l, n

    do n = 1, nv
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, ep
          do j = 1, ep
          do i = 1, ep
            pv(1,i,j,k,l,n) = real(ve(i,j,k,e,1,n), VTK_FLOAT64)
            pv(2,i,j,k,l,n) = real(ve(i,j,k,e,2,n), VTK_FLOAT64)
            pv(3,i,j,k,l,n) = real(ve(i,j,k,e,3,n), VTK_FLOAT64)
          end do
          end do
          end do
        end if
      end do
    end do

  end subroutine BuildLinearVectorData

  !-----------------------------------------------------------------------------
  !> Initializes the interpolation to quadratic cells.

  subroutine BuildInterpolationOperator(po, iop)
    integer,                intent(in)  :: po       !< polynomial order
    real(RNP), allocatable, intent(out) :: iop(:,:) !< 1D interpolation operator

    real(RNP) :: xc(0:po)   ! collocation points
    real(RNP) :: xi(2*po+1) ! interpolation points

    integer :: j, k, ni

    ni = size(xi)
    xc = LobattoPoints(po)
    xi(1::2) = xc
    xi(2::2) = HALF * (xc(0:po-1) + xc(1:po))

    allocate(iop(ni, 0:po))

    do k = 0, po
    do j = 1, ni
      iop(j,k) = LobattoPolynomial(k, xc, xi(j))
    end do
    end do

  end subroutine BuildInterpolationOperator

  !-----------------------------------------------------------------------------
  !> Interpolates element points to quadratic grid cells

  subroutine BuildQuadraticPointCoords(ep, ni, ne, nx, mask, iop, xe, points)
    integer,   intent(in)  :: ep         !< element points per direction
    integer,   intent(in)  :: ni         !< interpolated points per direction
    integer,   intent(in)  :: ne         !< num elements
    integer,   intent(in)  :: nx         !< number of exported elements
    logical,   intent(in)  :: mask(ne)   !< element mask
    real(RNP), intent(in)  :: iop(ni,ep) !< interpolation operator
    real(RNP), intent(in)  :: xe(ep,ep,ep,ne,3) !< mesh points
    real(VTK_FLOAT64), intent(out) :: points(3,ni,ni,ni,nx) !< VTK mesh points

    real(RNP), allocatable, save :: w(:,:,:,:)
    integer :: d, e, i, j, k, l

    !$omp master
    allocate(w(ni,ni,ni,ne))
    !$omp end master
    !$omp barrier

    do d = 1, 3
      call TPO_AAA(iop, xe(:,:,:,:,d), w)
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, ni
          do j = 1, ni
          do i = 1, ni
            points(d,i,j,k,l) = real(w(i,j,k,e), VTK_FLOAT64)
          end do
          end do
          end do
        end if
      end do
    end do

    !$omp master
    deallocate(w)
    !$omp end master

  end subroutine BuildQuadraticPointCoords

  !-----------------------------------------------------------------------------
  !> Connectivity of triquadratic hexahedral cells.

  subroutine BuildQuadraticCells(po, nx, cells)
    integer, intent(in) :: po  !< order of elements
    integer, intent(in) :: nx  !< number of exported elements
    integer(VTK_INT32), allocatable, intent(out) :: cells(:,:) !< VTK cells

    integer ::  c, i, j, k, l, m, n, o

    allocate(cells(27, nx * po**3))

    m = 2*po      ! grid intervals within one element
    n = (m+1)**3  ! grid points per element
    o = -1        ! offset of point IDs
    c =  1        ! cell counter

    do l = 1, nx
      do k = 1, 2*po, 2
      do j = 1, 2*po, 2
      do i = 1, 2*po, 2
        cells( 1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, m, m)
        cells( 2,c) = o + LexicalVertexIndex(i+1, j-1, k-1, m, m)
        cells( 3,c) = o + LexicalVertexIndex(i+1, j+1, k-1, m, m)
        cells( 4,c) = o + LexicalVertexIndex(i-1, j+1, k-1, m, m)
        cells( 5,c) = o + LexicalVertexIndex(i-1, j-1, k+1, m, m)
        cells( 6,c) = o + LexicalVertexIndex(i+1, j-1, k+1, m, m)
        cells( 7,c) = o + LexicalVertexIndex(i+1, j+1, k+1, m, m)
        cells( 8,c) = o + LexicalVertexIndex(i-1, j+1, k+1, m, m)
        cells( 9,c) = o + LexicalVertexIndex(i+0, j-1, k-1, m, m)
        cells(10,c) = o + LexicalVertexIndex(i+1, j+0, k-1, m, m)
        cells(11,c) = o + LexicalVertexIndex(i+0, j+1, k-1, m, m)
        cells(12,c) = o + LexicalVertexIndex(i-1, j+0, k-1, m, m)
        cells(13,c) = o + LexicalVertexIndex(i+0, j-1, k+1, m, m)
        cells(14,c) = o + LexicalVertexIndex(i+1, j+0, k+1, m, m)
        cells(15,c) = o + LexicalVertexIndex(i+0, j+1, k+1, m, m)
        cells(16,c) = o + LexicalVertexIndex(i-1, j+0, k+1, m, m)
        cells(17,c) = o + LexicalVertexIndex(i-1, j-1, k+0, m, m)
        cells(18,c) = o + LexicalVertexIndex(i+1, j-1, k+0, m, m)
        cells(19,c) = o + LexicalVertexIndex(i+1, j+1, k+0, m, m)
        cells(20,c) = o + LexicalVertexIndex(i-1, j+1, k+0, m, m)
        cells(21,c) = o + LexicalVertexIndex(i-1, j+0, k+0, m, m)
        cells(22,c) = o + LexicalVertexIndex(i+1, j+0, k+0, m, m)
        cells(23,c) = o + LexicalVertexIndex(i+0, j-1, k+0, m, m)
        cells(24,c) = o + LexicalVertexIndex(i+0, j+1, k+0, m, m)
        cells(25,c) = o + LexicalVertexIndex(i+0, j+0, k-1, m, m)
        cells(26,c) = o + LexicalVertexIndex(i+0, j+0, k+1, m, m)
        cells(27,c) = o + LexicalVertexIndex(i+0, j+0, k+0, m, m)
        c = c + 1
      end do
      end do
      end do
      o = o + n
    end do

  end subroutine BuildQuadraticCells

  !-----------------------------------------------------------------------------
  !> Interpolates scalar element variables to quadratic grid cells

  subroutine BuildQuadraticScalarData(ep, ni, ne, nx, ns, mask, iop, se, ps)
    integer,  intent(in)  :: ep         !< element points per direction
    integer,  intent(in)  :: ni         !< interpolated points per direction
    integer,  intent(in)  :: ne         !< num elements
    integer,  intent(in)  :: nx         !< number of exported elements
    integer,  intent(in)  :: ns         !< number of scalars
    logical,  intent(in)  :: mask(ne)   !< element mask
    real(RNP),intent(in)  :: iop(ni,ep) !< interpolation operator
    real(RNP),intent(in)  :: se(ep,ep,ep,ne,ns) !< SEM scalars
    real(VTK_FLOAT64), intent(out) :: ps(ni,ni,ni,nx,ns) !< VTK scalars

    real(RNP), allocatable, save :: w(:,:,:,:)
    integer :: e, i, j, k, l, n

    !$omp master
    allocate(w(ni,ni,ni,ne))
    !$omp end master
    !$omp barrier

    do n = 1, ns
      call TPO_AAA(iop, se(:,:,:,:,n), w)
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, ni
          do j = 1, ni
          do i = 1, ni
            ps(i,j,k,l,n) = real(w(i,j,k,e), VTK_FLOAT64)
          end do
          end do
          end do
        end if
      end do
    end do

    !$omp master
    deallocate(w)
    !$omp end master

  end subroutine BuildQuadraticScalarData

  !-----------------------------------------------------------------------------
  !> Interpolates vector element variables to quadratic grid cells

  subroutine BuildQuadraticVectorData(ep, ni, ne, nx, nv, mask, iop, ve, pv)
    integer,   intent(in)  :: ep         !< element points per direction
    integer,   intent(in)  :: ni         !< interpolated points per direction
    integer,   intent(in)  :: ne         !< num elements
    integer,   intent(in)  :: nx         !< number of exported elements
    integer,   intent(in)  :: nv         !< number of vectors
    logical,   intent(in)  :: mask(ne)   !< element mask
    real(RNP), intent(in)  :: iop(ni,ep) !< interpolation operator
    real(RNP), intent(in)  :: ve(ep,ep,ep,ne,3,nv) !< SEM scalars
    real(VTK_FLOAT64), intent(out) :: pv(3,ni,ni,ni,nx,nv) !< VTK scalars

    real(RNP), allocatable, save :: w(:,:,:,:)
    integer :: d, e, i, j, k, l, n

    !$omp master
    allocate(w(ni,ni,ni,ne))
    !$omp end master
    !$omp barrier

    do n = 1, nv
    do d = 1, 3
      call TPO_AAA(iop, ve(:,:,:,:,d,n), w)
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, ni
          do j = 1, ni
          do i = 1, ni
            pv(d,i,j,k,l,n) = real(w(i,j,k,e), VTK_FLOAT64)
          end do
          end do
          end do
        end if
      end do
    end do
    end do

    !$omp barrier
    !$omp master
    deallocate(w)
    !$omp end master

  end subroutine BuildQuadraticVectorData

  !=============================================================================

end module VTK__Export_Mesh_Data__3D
