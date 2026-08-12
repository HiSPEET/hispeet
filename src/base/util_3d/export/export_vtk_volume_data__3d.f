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
!> date:     2014/11/27, revised 2016-2021
!===============================================================================

module Export_VTK_Volume_Data__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Gauss_Jacobi,    only: LobattoPoints, LobattoPolynomial
  use C_Binding
  use VTK_Binding
  use TPO__AAA__3D
  use Structured_Mesh_Indexing__3D
  implicit none
  private

  public :: ExportVTK_VolumeData

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

  subroutine ExportVTK_VolumeData(x, s, sname, v, vname, file, part, n_parts, &
                                  subdiv, mask)

    ! mesh element collocation points
    real(RNP), intent(in) :: x(0:,0:,0:,:,:) !< mesh points [0:po,0:po,0:po,ne,3]

    ! scalar variables (optional)
    real(RNP), intent(in) :: s(0:,0:,0:,:,:) !< scalar [0:po,0:po,0:po,ne,ns]
    character(len=*), intent(in) :: sname(:) !< scalar names [ns]
    optional :: s, sname

    ! vector variables (optional)
    real(RNP), intent(in) :: v(0:,0:,0:,:,:,:) !< vectors [0:po,0:po,0:po,ne,3,nv]
    character(len=*), intent(in) :: vname(:)   !< vector names [nv]
    optional :: v, vname

    ! output control
    character(len=*),  intent(in) :: file    !< VTK output file
    integer, optional, intent(in) :: part    !< partition (piece)
    integer, optional, intent(in) :: n_parts !< number of partitions (pieces)
    logical, optional, intent(in) :: subdiv  !< T: quadratic subdivision [auto]
    logical, optional, intent(in) :: mask(:) !< skip elements with mask(e) = F

    ! VTK data .................................................................

    integer(C_INT) :: vtk         ! writer handle
    integer(C_INT) :: cell_type   ! cell type
    integer(C_INT) :: success     ! success flag

    integer(C_VTK_ID), allocatable :: cell(:,:)
    real(C_DOUBLE),    allocatable :: xg(:,:)
    real(C_DOUBLE),    allocatable :: sg(:,:)
    real(C_DOUBLE),    allocatable :: vg(:,:,:)
    logical,           allocatable :: mask_(:)

    ! auxiliary variables ......................................................

    real(RNP), allocatable :: iop(:,:) ! interpolation operator
    integer :: po, ne, ng, ni, nl, np, ns, nv
    integer :: interpolation_order, k
    character(len=80) :: tag

    ! initialization ...........................................................

    if (present(part)) then
      if (part < 0) then
        return ! skip empty partition
      else
        write(tag, fmt='(A2,I0)') '_p', part
      end if
    else
      tag = ''
    end if

    po = ubound(x,1)
    ne = ubound(x,4)

    if (present(mask)) then
      allocate(mask_(ne), source = mask)
    else
      allocate(mask_(ne), source = .true.)
    end if

    ! number of exported elements
    nl = count(mask_)

    ! number of scalars
    if (present(s)) then
      ns = size(s,5)
    else
      ns = 0
    end if

    ! number of vectors
    if (present(v)) then
      nv = size(v,6)
    else
      nv = 0
    end if

    ! automatic selection of interpolation order
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

    ! number of element points per direction
    np = po + 1
    if (interpolation_order == 2) then
      ! number of interpolated element points with quadratic interpolation
      ni = 2*po + 1
    else
      ni = np
    end if

    ! set up VTK file ..........................................................

    call VTK_XMLWriter_New(vtk)
    call VTK_XMLWriter_SetDataObjectType(vtk, VTK_UNSTRUCTURED_GRID)
    call VTK_XMLWriter_SetDataModeType(vtk, VTK_APPENDED)
    call VTK_XMLWriter_SetFileName(vtk, trim(file)//trim(tag)//'.vtu')

    ! grid points ..............................................................

    ng = nl * (interpolation_order*po + 1)**3

    allocate(xg(3,ng))

    select case(interpolation_order)
    case(1)
      call BuildLinearPointCoords(np, ne, nl, mask_, x, xg)
    case(2)
      call BuildInterpolationOperator(po, iop)
      call BuildQuadraticPointCoords(np, ni, ne, nl, mask_, iop, x, xg)
    end select

    call VTK_XMLWriter_SetPoints(vtk, xg)

    ! grid cells ...............................................................

    select case(interpolation_order)
    case(1)
      cell_type = VTK_HEXAHEDRON
      call BuildLinearCells(po, nl, cell)
    case(2)
      cell_type = VTK_TRIQUADRATIC_HEXAHEDRON
      call BuildQuadraticCells(po, nl, cell)
    end select

    call VTK_XMLWriter_SetCellsWithType(vtk, cell_type, cell)

    ! scalars ..................................................................

    if (present(s) .and. present(sname)) then
      allocate(sg(ng,ns))

      select case(interpolation_order)
      case(1)
        call BuildLinearScalarData(np, ne, nl, ns, mask_, s, sg)
      case(2)
        call BuildQuadraticScalarData(np, ni, ne, nl, ns, mask_, iop, s, sg)
      end select

      do k = 1, ns
        call VTK_XMLWriter_SetPointData(vtk, sname(k), sg(:,k), '')
      end do

    end if

    ! vectors ..................................................................

    if (present(v) .and. present(vname)) then
      allocate(vg(3,ng,nv))

      select case(interpolation_order)
      case(1)
        call BuildLinearVectorData(np, ne, nl, nv, mask_, v, vg)
      case(2)
        call BuildQuadraticVectorData(np, ni, ne, nl, nv, mask_, iop, v, vg)
      end select

      do k = 1, nv
        call VTK_XMLWriter_SetPointData(vtk, vname(k), vg(:,:,k), '')
      end do

    end if

    ! write ...................................................................

    call VTK_XMLWriter_Write(vtk, success)
    call VTK_XMLWriter_Delete(vtk)

    ! PVTU file ................................................................

    if (present(part) .and. present(n_parts)) then
      if (part == 0) then
        call Write_PVTU_File
      end if
    end if

  contains

    !---------------------------------------------------------------------------
    !> Writes the PVTU file for a parallel multi-piece data set

    subroutine Write_PVTU_File

      integer :: pvtu

      open(newunit=pvtu, file=trim(file)//'.pvtu')

      ! header
      write(pvtu,'(1A)') '<?xml version="1.0"?>'
      write(pvtu,'(3A)') '<VTKFile type="PUnstructuredGrid" version="0.1" ', &
                         'byte_order="LittleEndian" ',                       &
                         'compressor="vtkZLibDataCompressor">'

      write(pvtu,'(2X,A)') '<PUnstructuredGrid GhostLevel="0">'

      ! start point data section
      write(pvtu,'(4X,A)') '<PPointData>'

      ! scalars
      do k = 1, ns
        write(pvtu,'(6X,4A)') '<PDataArray type="Float64" ', &
                              'Name="', trim(sname(k)), '"/>'
      end do

      ! vectors
      do k = 1, nv
        write(pvtu,'(6X,5A)') '<PDataArray type="Float64" ', &
                              'Name="', trim(vname(k)),'" ', &
                              'NumberOfComponents="3"/>'
      end do

      ! close point data section
      write(pvtu,'(4X,A)') '</PPointData>'

      ! points section
      write(pvtu,'(4X,A)') '<PPoints>'
      write(pvtu,'(6X,A)') '<PDataArray type="Float64" NumberOfComponents="3"/>'
      write(pvtu,'(4X,A)') '</PPoints>'

      ! piece sources
      do k = 0, n_parts-1
        write(pvtu,'(4X,3A,I0,A)') '<Piece Source="',trim(file),'_p',k,'.vtu"/>'
      end do

      ! trailer
      write(pvtu,'(2X,1A)') '</PUnstructuredGrid>'
      write(pvtu,'(A)') '</VTKFile>'

      ! close file
      close(pvtu)

    end subroutine Write_PVTU_File

  end subroutine ExportVTK_VolumeData

  !-----------------------------------------------------------------------------
  !> Maps element points to linear grid cells

  subroutine BuildLinearPointCoords(np, ne, nl, mask, xe, xg)
    integer,        intent(in)  :: np       !< element points per direction
    integer,        intent(in)  :: ne       !< num elements
    integer,        intent(in)  :: nl       !< num exported elements
    logical,        intent(in)  :: mask(ne) !< element mask
    real(RNP),      intent(in)  :: xe(np,np,np,ne,3) !< mesh points
    real(C_DOUBLE), intent(out) :: xg(3,np,np,np,nl) !< VTK grid points

    integer :: e, i, j, k, l

    l = 0
    do e = 1, ne
      if (mask(e)) then
        l = l + 1
        do k = 1, np
        do j = 1, np
        do i = 1, np
          xg(1,i,j,k,l) = real(xe(i,j,k,e,1), C_DOUBLE)
          xg(2,i,j,k,l) = real(xe(i,j,k,e,2), C_DOUBLE)
          xg(3,i,j,k,l) = real(xe(i,j,k,e,3), C_DOUBLE)
        end do
        end do
        end do
      end if
    end do

  end subroutine BuildLinearPointCoords

  !-----------------------------------------------------------------------------
  !> Connectivity of trilinear hexahedral cells.

  subroutine BuildLinearCells(po, nl, cell)
    integer, intent(in) :: po  !< order of elements
    integer, intent(in) :: nl  !< number of exported elements
    integer(C_VTK_ID), allocatable, intent(out) :: cell(:,:) !< grid cells

    integer :: c, i, j, k, l, n, o

    allocate(cell(0:8, nl * po**3))

    cell(0,:) = 8    ! grid points per cell
    n = (po + 1)**3  ! grid points per element
    o = -1           ! offset of point IDs
    c =  1           ! cell counter

    do l = 1, nl
      do k = 1, po
      do j = 1, po
      do i = 1, po
        cell(1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, po, po)
        cell(2,c) = o + LexicalVertexIndex(i  , j-1, k-1, po, po)
        cell(3,c) = o + LexicalVertexIndex(i  , j  , k-1, po, po)
        cell(4,c) = o + LexicalVertexIndex(i-1, j  , k-1, po, po)
        cell(5,c) = o + LexicalVertexIndex(i-1, j-1, k  , po, po)
        cell(6,c) = o + LexicalVertexIndex(i  , j-1, k  , po, po)
        cell(7,c) = o + LexicalVertexIndex(i  , j  , k  , po, po)
        cell(8,c) = o + LexicalVertexIndex(i-1, j  , k  , po, po)
        c = c + 1
      end do
      end do
      end do
      o = o + n
    end do

  end subroutine BuildLinearCells

  !-----------------------------------------------------------------------------
  !> Maps scalar element variables to linear grid cells

  subroutine BuildLinearScalarData(np, ne, nl, ns, mask, se, sg)
    integer,        intent(in)  :: np       !< element points per direction
    integer,        intent(in)  :: ne       !< num elements
    integer,        intent(in)  :: nl       !< num exported elements
    integer,        intent(in)  :: ns       !< num scalars
    logical,        intent(in)  :: mask(ne) !< element mask
    real(RNP),      intent(in)  :: se(np,np,np,ne,ns) !< SEM scalars
    real(C_DOUBLE), intent(out) :: sg(np,np,np,nl,ns) !< VTK scalars

    integer :: e, i, j, k, l, n

    do n = 1, ns
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, np
          do j = 1, np
          do i = 1, np
            sg(i,j,k,l,n) = real(se(i,j,k,e,n), C_DOUBLE)
          end do
          end do
          end do
        end if
      end do
    end do

  end subroutine BuildLinearScalarData

  !-----------------------------------------------------------------------------
  !> Maps vector element variables to linear grid cells

  subroutine BuildLinearVectorData(np, ne, nl, nv, mask, ve, vg)
    integer,        intent(in)  :: np       !< element points per direction
    integer,        intent(in)  :: ne       !< num elements
    integer,        intent(in)  :: nl       !< num exported elements
    integer,        intent(in)  :: nv       !< num vectors
    logical,        intent(in)  :: mask(ne) !< element mask
    real(RNP),      intent(in)  :: ve(np,np,np,ne,3,nv) !< SEM vectors
    real(C_DOUBLE), intent(out) :: vg(3,np,np,np,nl,nv) !< VTK vectors

    integer :: e, i, j, k, l, n

    do n = 1, nv
      l = 0
      do e = 1, ne
        if (mask(e)) then
          l = l + 1
          do k = 1, np
          do j = 1, np
          do i = 1, np
            vg(1,i,j,k,l,n) = real(ve(i,j,k,e,1,n), C_DOUBLE)
            vg(2,i,j,k,l,n) = real(ve(i,j,k,e,2,n), C_DOUBLE)
            vg(3,i,j,k,l,n) = real(ve(i,j,k,e,3,n), C_DOUBLE)
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

  subroutine BuildQuadraticPointCoords(np, ni, ne, nl, mask, iop, xe, xg)
    integer,        intent(in)  :: np         !< element points per direction
    integer,        intent(in)  :: ni         !< interpolated points per direct.
    integer,        intent(in)  :: ne         !< num elements
    integer,        intent(in)  :: nl         !< number of exported elements
    logical,        intent(in)  :: mask(ne)   !< element mask
    real(RNP),      intent(in)  :: iop(ni,np) !< interpolation operator
    real(RNP),      intent(in)  :: xe(np,np,np,ne,3) !< mesh points
    real(C_DOUBLE), intent(out) :: xg(3,ni,ni,ni,nl) !< VTK grid points

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
            xg(d,i,j,k,l) = real(w(i,j,k,e), C_DOUBLE)
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

  subroutine BuildQuadraticCells(po, nl, cell)
    integer, intent(in) :: po  !< order of elements
    integer, intent(in) :: nl  !< number of exported elements
    integer(C_VTK_ID), allocatable, intent(out) :: cell(:,:) !< grid cells

    integer ::  c, i, j, k, l, m, n, o

    allocate(cell(0:27, nl * po**3))

    cell(0,:) = 27   ! grid points per cell
    m = 2*po         ! grid intervals within one element
    n = (m+1)**3     ! grid points per element
    o = -1           ! offset of point IDs
    c =  1           ! cell counter

    do l = 1, nl
      do k = 1, 2*po, 2
      do j = 1, 2*po, 2
      do i = 1, 2*po, 2
        cell( 1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, m, m)
        cell( 2,c) = o + LexicalVertexIndex(i+1, j-1, k-1, m, m)
        cell( 3,c) = o + LexicalVertexIndex(i+1, j+1, k-1, m, m)
        cell( 4,c) = o + LexicalVertexIndex(i-1, j+1, k-1, m, m)
        cell( 5,c) = o + LexicalVertexIndex(i-1, j-1, k+1, m, m)
        cell( 6,c) = o + LexicalVertexIndex(i+1, j-1, k+1, m, m)
        cell( 7,c) = o + LexicalVertexIndex(i+1, j+1, k+1, m, m)
        cell( 8,c) = o + LexicalVertexIndex(i-1, j+1, k+1, m, m)
        cell( 9,c) = o + LexicalVertexIndex(i+0, j-1, k-1, m, m)
        cell(10,c) = o + LexicalVertexIndex(i+1, j+0, k-1, m, m)
        cell(11,c) = o + LexicalVertexIndex(i+0, j+1, k-1, m, m)
        cell(12,c) = o + LexicalVertexIndex(i-1, j+0, k-1, m, m)
        cell(13,c) = o + LexicalVertexIndex(i+0, j-1, k+1, m, m)
        cell(14,c) = o + LexicalVertexIndex(i+1, j+0, k+1, m, m)
        cell(15,c) = o + LexicalVertexIndex(i+0, j+1, k+1, m, m)
        cell(16,c) = o + LexicalVertexIndex(i-1, j+0, k+1, m, m)
        cell(17,c) = o + LexicalVertexIndex(i-1, j-1, k+0, m, m)
        cell(18,c) = o + LexicalVertexIndex(i+1, j-1, k+0, m, m)
        cell(19,c) = o + LexicalVertexIndex(i+1, j+1, k+0, m, m)
        cell(20,c) = o + LexicalVertexIndex(i-1, j+1, k+0, m, m)
        cell(21,c) = o + LexicalVertexIndex(i-1, j+0, k+0, m, m)
        cell(22,c) = o + LexicalVertexIndex(i+1, j+0, k+0, m, m)
        cell(23,c) = o + LexicalVertexIndex(i+0, j-1, k+0, m, m)
        cell(24,c) = o + LexicalVertexIndex(i+0, j+1, k+0, m, m)
        cell(25,c) = o + LexicalVertexIndex(i+0, j+0, k-1, m, m)
        cell(26,c) = o + LexicalVertexIndex(i+0, j+0, k+1, m, m)
        cell(27,c) = o + LexicalVertexIndex(i+0, j+0, k+0, m, m)
        c = c + 1
      end do
      end do
      end do
      o = o + n
    end do

  end subroutine BuildQuadraticCells

  !-----------------------------------------------------------------------------
  !> Interpolates scalar element variables to quadratic grid cells

  subroutine BuildQuadraticScalarData(np, ni, ne, nl, ns, mask, iop, se, sg)
    integer,        intent(in)  :: np         !< element points per direction
    integer,        intent(in)  :: ni         !< interpolated points per direct.
    integer,        intent(in)  :: ne         !< num elements
    integer,        intent(in)  :: nl         !< number of exported elements
    integer,        intent(in)  :: ns         !< number of scalars
    logical,        intent(in)  :: mask(ne)   !< element mask
    real(RNP),      intent(in)  :: iop(ni,np) !< interpolation operator
    real(RNP),      intent(in)  :: se(np,np,np,ne,ns) !< SEM scalars
    real(C_DOUBLE), intent(out) :: sg(ni,ni,ni,nl,ns) !< VTK scalars

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
            sg(i,j,k,l,n) = real(w(i,j,k,e), C_DOUBLE)
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

  subroutine BuildQuadraticVectorData(np, ni, ne, nl, nv, mask, iop, ve, vg)
    integer,        intent(in)  :: np         !< element points per direction
    integer,        intent(in)  :: ni         !< interpolated points per direct.
    integer,        intent(in)  :: ne         !< num elements
    integer,        intent(in)  :: nl         !< number of exported elements
    integer,        intent(in)  :: nv         !< number of vectors
    logical,        intent(in)  :: mask(ne)   !< element mask
    real(RNP),      intent(in)  :: iop(ni,np) !< interpolation operator
    real(RNP),      intent(in)  :: ve(np,np,np,ne,3,nv) !< SEM scalars
    real(C_DOUBLE), intent(out) :: vg(3,ni,ni,ni,nl,nv) !< VTK scalars

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
            vg(d,i,j,k,l,n) = real(w(i,j,k,e), C_DOUBLE)
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

end module Export_VTK_Volume_Data__3D
