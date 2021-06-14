!> summary:  Export mesh data into VTK XML file
!> author:   Joerg Stiller
!> date:     2014/11/27, revised 2016/12/02
!> license:  Institute of Fluid Mechanics TU Dresden, 01062 Dresden, Germany
!>
!>### Export mesh data into VTK XML file
!===============================================================================

module Export_Volume_Data_To_VTK
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Gauss_Jacobi,    only: LobattoPoints, LobattoPolynomial
  use C_Binding
  use VTK_Binding
  use TPO__AAA__3D
  use Structured_Mesh  ! adapt vertex numbering when switching to Mesh_3d__Structured_Indexing
  implicit none
  private

  public :: ExportVolumeDataToVTK

contains

!-------------------------------------------------------------------------------
!> Export elements and variables of mesh partition into VTK XML file.
!>
!> In case of partitioned (multi-piece) data, the partition ID, part, and the
!> number of partitions, n_part, must be given. The partitions are numbered
!> from 0 to n_part-1. Partition 0 writes the PVTU (parallel unstructured mesh)
!> file.

subroutine ExportVolumeDataToVTK(po, ne, ns, nv, x, s, sname, v, vname, file, &
    part, n_part)

  ! dimensions
  integer, intent(in) :: po  !< polynomial order
  integer, intent(in) :: ne  !< number of elements
  integer, intent(in) :: ns  !< number of scalars
  integer, intent(in) :: nv  !< number of vectors

  ! mesh element collocation points
  real(RNP), intent(in) :: x(0:po,0:po,0:po,ne,3)           !< mesh points

  ! scalars (optional)
  real(RNP),        intent(in) :: s(0:po,0:po,0:po,ne,ns)   !< scalar variables
  character(len=*), intent(in) :: sname(ns)                 !< scalar names
  optional :: s, sname

  ! vectors (optional)
  real(RNP),        intent(in) :: v(0:po,0:po,0:po,ne,3,nv) !< vector variables
  character(len=*), intent(in) :: vname(nv)                 !< vector names
  optional :: v, vname

  ! output control
  character(len=*),  intent(in) :: file   !< VTK output file
  integer, optional, intent(in) :: part   !< partition (piece)
  integer, optional, intent(in) :: n_part !< number of partitions (pieces)

  ! VTK data ...................................................................

  integer(C_INT) :: vtk         ! writer handle
  integer(C_INT) :: cell_type   ! cell type
  integer(C_INT) :: success     ! success flag

  integer(C_VTK_ID), allocatable :: cell(:,:)
  real(C_DOUBLE),    allocatable :: xg(:,:)
  real(C_DOUBLE),    allocatable :: sg(:,:)
  real(C_DOUBLE),    allocatable :: vg(:,:,:)

  ! auxiliary variables ........................................................

  character(len=80) :: tag
  integer :: interpolation_order, np, k
  real(RNP), allocatable :: iop(:,:) ! interpolation operator

  ! initialization .............................................................

  ! automatic selection of interpolation order
  if (po == 1) then
    interpolation_order = 1
  else if (po > 1) then
    interpolation_order = 2
  else
    return
  end if

  if (present(part)) then
    write(tag, fmt='(A2,I0)') '_p', part
  else
    tag = ''
  end if

  ! set up VTK file ............................................................

  call VTK_XMLWriter_New(vtk)
  call VTK_XMLWriter_SetDataObjectType(vtk, VTK_UNSTRUCTURED_GRID)
  call VTK_XMLWriter_SetDataModeType(vtk, VTK_APPENDED)
  call VTK_XMLWriter_SetFileName(vtk, trim(file)//trim(tag)//'.vtu')

  ! grid points ................................................................

  np = ne * (interpolation_order*po + 1)**3

  allocate(xg(3,np))

  select case(interpolation_order)
  case(1)
    call BuildLinearPointCoords(np, x, xg)
  case(2)
    call BuildInterpolationOperator(po, iop)
    call BuildQuadraticPointCoords(iop, x, xg)
  end select

  call VTK_XMLWriter_SetPoints(vtk, xg)

  ! grid cells .................................................................

  select case(interpolation_order)
  case(1)
    cell_type = VTK_HEXAHEDRON
    call BuildLinearCells(po, ne, cell)
  case(2)
    cell_type = VTK_TRIQUADRATIC_HEXAHEDRON
    call BuildQuadraticCells(po, ne, cell)
  end select

  call VTK_XMLWriter_SetCellsWithType(vtk, cell_type, cell)

  ! scalars ....................................................................

  if (present(s) .and. present(sname)) then
    allocate(sg(np,ns))

    select case(interpolation_order)
    case(1)
      call BuildLinearScalarData(s, sg)
    case(2)
      call BuildQuadraticScalarData(iop, s, sg)
    end select

    do k = 1, ns
      call VTK_XMLWriter_SetPointData(vtk, sname(k), sg(:,k), '')
    end do

  end if

  ! vectors ....................................................................

  if (present(v) .and. present(vname)) then
    allocate(vg(3,np,nv))

    select case(interpolation_order)
    case(1)
      call BuildLinearVectorData(np, nv, v, vg)
    case(2)
      call BuildQuadraticVectorData(iop, v, vg)
    end select

    do k = 1, nv
      call VTK_XMLWriter_SetPointData(vtk, vname(k), vg(:,:,k), '')
    end do

  end if

  ! write ......................................................................

  call VTK_XMLWriter_Write(vtk, success)
  call VTK_XMLWriter_Delete(vtk)

  ! PVTU file ..................................................................

  if (present(part) .and. present(n_part)) then
    if (part == 0) then
      call Write_PVTU_File
    end if
  end if

contains

  !-----------------------------------------------------------------------------
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

    ! close poit data section
    write(pvtu,'(4X,A)') '</PPointData>'

    ! points section
    write(pvtu,'(4X,A)') '<PPoints>'
    write(pvtu,'(6X,A)') '<PDataArray type="Float64" NumberOfComponents="3"/>'
    write(pvtu,'(4X,A)') '</PPoints>'

    ! piece sources
    do k = 0, n_part-1
      write(pvtu,'(4X,3A,I0,A)') '<Piece Source="',trim(file),'_p',k,'.vtu"/>'
    end do

    ! trailer
    write(pvtu,'(2X,1A)') '</PUnstructuredGrid>'
    write(pvtu,'(A)') '</VTKFile>'

    ! close file
    close(pvtu)

  end subroutine Write_PVTU_File

end subroutine ExportVolumeDataToVTK

!-------------------------------------------------------------------------------
!> Maps element points to linear grid cells

subroutine BuildLinearPointCoords(np, xc, xg)
  integer,        intent(in)  :: np       !< number of mesh points
  real(RNP),      intent(in)  :: xc(np,3) !< mesh element points
  real(C_DOUBLE), intent(out) :: xg(3,np) !< VTK grid points

  xg = transpose(xc)

end subroutine BuildLinearPointCoords

!-------------------------------------------------------------------------------
!> Connectivity of trilinear hexahedral cells.

subroutine BuildLinearCells(po, ne, cell)
  integer, intent(in) :: po  !< order of elements
  integer, intent(in) :: ne  !< number of elements
  integer(C_VTK_ID), allocatable, intent(out) :: cell(:,:) !< grid cells

  integer :: c, i, j, k, l, n, o

  allocate(cell(0:8, ne * po**3))

  cell(0,:) = 8    ! grid points per cell
  n = (po + 1)**3  ! grid points per element
  o = -1           ! offset of point IDs
  c =  1           ! cell counter

  do l = 1, ne
    do k = 1, po
    do j = 1, po
    do i = 1, po
      cell(1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, po, po, 0)
      cell(2,c) = o + LexicalVertexIndex(i  , j-1, k-1, po, po, 0)
      cell(3,c) = o + LexicalVertexIndex(i  , j  , k-1, po, po, 0)
      cell(4,c) = o + LexicalVertexIndex(i-1, j  , k-1, po, po, 0)
      cell(5,c) = o + LexicalVertexIndex(i-1, j-1, k  , po, po, 0)
      cell(6,c) = o + LexicalVertexIndex(i  , j-1, k  , po, po, 0)
      cell(7,c) = o + LexicalVertexIndex(i  , j  , k  , po, po, 0)
      cell(8,c) = o + LexicalVertexIndex(i-1, j  , k  , po, po, 0)
      c = c + 1
    end do
    end do
    end do
    o = o + n
  end do

end subroutine BuildLinearCells

!-------------------------------------------------------------------------------
!> Maps scalar element variables to linear grid cells

subroutine BuildLinearScalarData(sc, sg)
  real(RNP),      intent(in)  :: sc(:,:,:,:,:) !< scalars at collocation points
  real(C_DOUBLE), intent(out) :: sg(:,:)       !< scalars at VTK grid points

  sg = reshape(sc, shape(sg))

end subroutine BuildLinearScalarData

!-------------------------------------------------------------------------------
!> Maps vector element variables to linear grid cells

subroutine BuildLinearVectorData(np, nv, vc, vg)
  integer,        intent(in)  :: np          !< number of mesh points
  integer,        intent(in)  :: nv          !< number of vectors
  real(RNP),      intent(in)  :: vc(np,3,nv) !< vectors at collocation pts.
  real(C_DOUBLE), intent(out) :: vg(3,np,nv) !< vectors at VTK grid points

  integer :: i

  do i = 1, nv
    vg(:,:,i) = transpose(vc(:,:,i))
  end do

end subroutine BuildLinearVectorData

!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
!> Interpolates element points to quadratic grid cells

subroutine BuildQuadraticPointCoords(iop, xc, xg)
  real(RNP),      intent(in)  :: iop(:,:)      !< interpolation operator
  real(RNP),      intent(in)  :: xc(:,:,:,:,:) !< mesh element points
  real(C_DOUBLE), intent(out) :: xg(:,:)       !< VTK grid points

  real(RNP), allocatable, save :: w(:,:,:,:)
  integer :: ne, ng, np
  integer :: i

  ng = size(iop,1)
  ne = size(xc,4)
  np = size(xg,2)

  !$omp single
  allocate(w(ng,ng,ng,ne))
  !$omp end single

  do i = 1, 3
    call TPO_AAA(iop, xc(:,:,:,:,i), w)
    xg(i,:) = reshape(w, [np])
  end do

  !$omp barrier
  !$omp master
  deallocate(w)
  !$omp end master

end subroutine BuildQuadraticPointCoords

!-------------------------------------------------------------------------------
!> Connectivity of triquadratic hexahedral cells.

subroutine BuildQuadraticCells(po, ne, cell)
  integer, intent(in) :: po  !< order of elements
  integer, intent(in) :: ne  !< number of elements
  integer(C_VTK_ID), allocatable, intent(out) :: cell(:,:) !< grid cells

  integer ::  c, i, j, k, l, m, n, o

  allocate(cell(0:27, ne * po**3))

  cell(0,:) = 27   ! grid points per cell
  m = 2*po         ! grid intervals within one element
  n = (m+1)**3     ! grid points per element
  o = -1           ! offset of point IDs
  c =  1           ! cell counter

  do l = 1, ne
    do k = 1, 2*po, 2
    do j = 1, 2*po, 2
    do i = 1, 2*po, 2
      cell( 1,c) = o + LexicalVertexIndex(i-1, j-1, k-1, m, m, 0)
      cell( 2,c) = o + LexicalVertexIndex(i+1, j-1, k-1, m, m, 0)
      cell( 3,c) = o + LexicalVertexIndex(i+1, j+1, k-1, m, m, 0)
      cell( 4,c) = o + LexicalVertexIndex(i-1, j+1, k-1, m, m, 0)
      cell( 5,c) = o + LexicalVertexIndex(i-1, j-1, k+1, m, m, 0)
      cell( 6,c) = o + LexicalVertexIndex(i+1, j-1, k+1, m, m, 0)
      cell( 7,c) = o + LexicalVertexIndex(i+1, j+1, k+1, m, m, 0)
      cell( 8,c) = o + LexicalVertexIndex(i-1, j+1, k+1, m, m, 0)
      cell( 9,c) = o + LexicalVertexIndex(i+0, j-1, k-1, m, m, 0)
      cell(10,c) = o + LexicalVertexIndex(i+1, j+0, k-1, m, m, 0)
      cell(11,c) = o + LexicalVertexIndex(i+0, j+1, k-1, m, m, 0)
      cell(12,c) = o + LexicalVertexIndex(i-1, j+0, k-1, m, m, 0)
      cell(13,c) = o + LexicalVertexIndex(i+0, j-1, k+1, m, m, 0)
      cell(14,c) = o + LexicalVertexIndex(i+1, j+0, k+1, m, m, 0)
      cell(15,c) = o + LexicalVertexIndex(i+0, j+1, k+1, m, m, 0)
      cell(16,c) = o + LexicalVertexIndex(i-1, j+0, k+1, m, m, 0)
      cell(17,c) = o + LexicalVertexIndex(i-1, j-1, k+0, m, m, 0)
      cell(18,c) = o + LexicalVertexIndex(i+1, j-1, k+0, m, m, 0)
      cell(19,c) = o + LexicalVertexIndex(i+1, j+1, k+0, m, m, 0)
      cell(20,c) = o + LexicalVertexIndex(i-1, j+1, k+0, m, m, 0)
      cell(21,c) = o + LexicalVertexIndex(i-1, j+0, k+0, m, m, 0)
      cell(22,c) = o + LexicalVertexIndex(i+1, j+0, k+0, m, m, 0)
      cell(23,c) = o + LexicalVertexIndex(i+0, j-1, k+0, m, m, 0)
      cell(24,c) = o + LexicalVertexIndex(i+0, j+1, k+0, m, m, 0)
      cell(25,c) = o + LexicalVertexIndex(i+0, j+0, k-1, m, m, 0)
      cell(26,c) = o + LexicalVertexIndex(i+0, j+0, k+1, m, m, 0)
      cell(27,c) = o + LexicalVertexIndex(i+0, j+0, k+0, m, m, 0)
      c = c + 1
    end do
    end do
    end do
    o = o + n
  end do

end subroutine BuildQuadraticCells

!-------------------------------------------------------------------------------
!> Interpolates scalar element variables to quadratic grid cells

subroutine BuildQuadraticScalarData(iop, sc, sg)
  real(RNP),      intent(in)  :: iop(:,:)      !< interpolation operator
  real(RNP),      intent(in)  :: sc(:,:,:,:,:) !< scalars at collocation points
  real(C_DOUBLE), intent(out) :: sg(:,:)       !< scalars at VTK grid points

  real(RNP), allocatable, save :: w(:,:,:,:)
  integer :: ne, ng, np, ns
  integer :: i

  ng = size(iop,1)
  ne = size(sc,4)
  ns = size(sc,5)
  np = size(sg,1)

  !$omp single
  allocate(w(ng,ng,ng,ne))
  !$omp end single

  do i = 1, ns
    call TPO_AAA(iop, sc(:,:,:,:,i), w)
    sg(:,i) = reshape(w, [np])
  end do

  !$omp barrier
  !$omp master
  deallocate(w)
  !$omp end master

end subroutine BuildQuadraticScalarData

!-------------------------------------------------------------------------------
!> Interpolates vector element variables to quadratic grid cells

subroutine BuildQuadraticVectorData(iop, vc, vg)
  real(RNP),      intent(in)  :: iop(:,:)        !< interpolation operator
  real(RNP),      intent(in)  :: vc(:,:,:,:,:,:) !< vectors at collocation pts.
  real(C_DOUBLE), intent(out) :: vg(:,:,:)       !< vectors at VTK grid points

  real(RNP), allocatable, save :: w(:,:,:,:)
  integer :: ne, ng, np, nv
  integer :: i, j

  ng = size(iop,1)
  ne = size(vc,4)
  nv = size(vc,6)
  np = size(vg,2)

  !$omp single
  allocate(w(ng,ng,ng,ne))
  !$omp end single

  do j = 1, nv
    do i = 1, 3
      call TPO_AAA(iop, vc(:,:,:,:,i,j), w)
      vg(i,:,j) = reshape(w, [np])
    end do
  end do

  !$omp barrier
  !$omp master
  deallocate(w)
  !$omp end master

end subroutine BuildQuadraticVectorData

!===============================================================================

end module Export_Volume_Data_To_VTK
