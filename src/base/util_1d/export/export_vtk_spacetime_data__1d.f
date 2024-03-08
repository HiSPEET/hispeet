!> summary:  Export 1D mesh space-time data into VTK XML file
!> author:   Erik Pfister
!> date:     2023/09/25
!> license:  Institute of Fluid Mechanics TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Export_VTK_Spacetime_Data__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Gauss_Jacobi,    only: LobattoPoints, LobattoPolynomial
  use C_Binding
  use Execution_Control
  use VTK_Binding
  use TPO__AAA__3D
  use Structured_Mesh_Indexing__3D
  implicit none
  private

  public :: Pack_SpacetimeData
  public :: ExportVTK_SpacetimeData

contains

  !-----------------------------------------------------------------------------
  !> Packs variables of the same spacetime grid in one scalar array.
  !>
  !> This routine takes variables 'u' with corresponding names 'uname' and
  !> extends the arrays 's' and 'sname'. This routine can be called multiple
  !> times before the routine ExportVTK_SpacetimeData is called.

  subroutine Pack_SpacetimeData(u, s, uname, sname)

    real(RNP), intent(in) :: u(:,:,:,:,:)
      !< variables to be appended to the array `s` (0:po,ne,nc,0:pt,nt)
    real(RNP), allocatable, intent(inout) :: s(:,:,:,:,:)
      !< array of scalars to be extended (0:po,ne,0:pt,nt,ns->ns_new)
    character(len=*), intent(in) :: uname(:)
      !< names of the variables appended (nc)
    character(len=20), allocatable, intent(inout) :: sname(:)
      !< names of scalars contained in `s`(ns->ns_new)

    ! auxiliary variables ......................................................

    real(RNP)        , allocatable :: temp_s(:,:,:,:,:)
    character(len=20), allocatable :: temp_sname(:)

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
      call Error( 'Pack_SpacetimeData'                                       &
                , 'Length of uname does not match number of components in u' &
                , 'Export_VTK_Spacetime_Data__1D'                            )
    endif

    ! update ns
    ns_new = ns + nc

    ! check if s is not already allocated, allocate it with the correct dimensions
    if (.not. allocated(s)) then
      allocate(s(ubound(u,1),ubound(u,2),ubound(u,4),ubound(u,5),ns_new))
      allocate(sname(ns_new))
    else
      ! reallocate s to accommodate the new variables
      allocate(temp_s(ubound(s,1),ubound(s,2),ubound(s,3),ubound(s,4),ns_new))
      allocate(temp_sname(ns_new))
      temp_s(:,:,:,:,1:ns) = s
      temp_sname(1:ns)     = sname
      call move_alloc(temp_s,     s    )
      call move_alloc(temp_sname, sname)
    endif

    ! append new scalar variables for s
    do i = 1, nc
      s(:,:,:,:,ns+i) = u(:,:,i,:,:)
      sname(ns+i)     = uname(i)
    enddo

  end subroutine Pack_SpacetimeData

  !-----------------------------------------------------------------------------
  !> Export elements and variables of mesh partition into VTK XML file.
  !>
  !> It is assumed that the mesh points are the element Lobatto nodes.
  !>
  !> This routine creates .vtu files of variables on a space-time grid
  !> for visualizations on ParaView. One axis represtens a spatial dimension
  !> and the other axis represents the evolution of that variable in time.

  subroutine ExportVTK_SpacetimeData(x, t, s, sname, file, part, n_parts, subdiv)
    real(RNP), intent(in) :: x(0:,:) !< mesh points in space (0:po,ne)
    real(RNP), intent(in) :: t(0:,:) !< time points (0:pt,1:nt)
    real(RNP), intent(in) :: s(0:,:,0:,:,:)   !< scalars (0:po,ne,0:pt,nt,ns)
    character(len=*),  intent(in) :: sname(:) !< scalar names (ns)
    character(len=*),  intent(in) :: file     !< VTK output file
    integer, optional, intent(in) :: part     !< partition (piece)
    integer, optional, intent(in) :: n_parts  !< number of partitions (pieces)
    logical, optional, intent(in) :: subdiv   !< T: quadratic subdivision [auto]
    optional :: s, sname

    ! VTK data .................................................................

    integer(C_INT) :: vtk         ! writer handle
    integer(C_INT) :: cell_type   ! cell type
    integer(C_INT) :: success     ! success flag

    integer(C_VTK_ID), allocatable :: cell(:,:)
    real(C_DOUBLE),    allocatable :: xg(:,:)
    real(C_DOUBLE),    allocatable :: sg(:,:)

    ! auxiliary variables ......................................................

    integer :: po, ne, nc, pt, nt, ns
    integer :: np, k
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
    ne = ubound(x,2)
    pt = ubound(t,1)
    nt = ubound(t,2)

    if (present(s)) then
      ns = size(s,5)
    else
      ns = 0
    end if

    ! set up VTK file ..........................................................

    ! this will map the spatial dimension against the time dimension

    call VTK_XMLWriter_New(vtk)
    call VTK_XMLWriter_SetDataObjectType(vtk, VTK_UNSTRUCTURED_GRID)
    call VTK_XMLWriter_SetDataModeType(vtk, VTK_APPENDED)
    call VTK_XMLWriter_SetFileName(vtk, trim(file)//trim(tag)//'.vtu')

    ! grid points ..............................................................

    ! for a 2D space-time grid, calculate the number of spatial points and
    ! multiply by the number of time steps
    np = ne * (po + 1) * nt * (pt + 1)

    allocate(xg(3,np))

    call BuildLinearSpacetimeCoords(np, po, ne, pt, nt, x, t, xg)

    call VTK_XMLWriter_SetPoints(vtk, xg)

    ! grid cells ...............................................................

    cell_type = VTK_QUAD
    call BuildLinearSpacetimeCells(po, pt, ne, nt, cell)

    call VTK_XMLWriter_SetCellsWithType(vtk, cell_type, cell)

    ! scalars ..................................................................

    if (present(s) .and. present(sname)) then
      allocate(sg(np,ns))

      call BuildLinearScalarData(np, ns, s, sg)

      do k = 1, ns
        call VTK_XMLWriter_SetPointData(vtk, sname(k), sg(:,k), '')
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

      ! close poit data section
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

  end subroutine ExportVTK_SpacetimeData

  !-----------------------------------------------------------------------------
  !> Maps spacetime element points to linear grid cells

  subroutine BuildLinearSpacetimeCoords(np, po, ne, pt, nt, xc, tc, xg)
    integer,        intent(in)  :: np          !< number of mesh points
    real(RNP),      intent(in)  :: xc(0:po,ne) !< space mesh element points
    real(RNP),      intent(in)  :: tc(0:pt,nt) !< time mesh element points
    real(C_DOUBLE), intent(out) :: xg(3,np)    !< VTK grid points

    integer :: i, j, k, l, idx, po, ne, pt, nt

    idx = 1
    do l = 1, nt
      do k = 0, pt
        do j = 1, ne
          do i = 0, po
            xg(1, idx) = xc(i, j)
            xg(2, idx) = tc(k, l)
            xg(3, idx) = 0
            idx = idx + 1
          end do
        end do
      end do
    end do

  end subroutine BuildLinearSpacetimeCoords

  !-----------------------------------------------------------------------------
  !> Connectivity of bilinear quadrilateral cells for spacetime grid.

  subroutine BuildLinearSpacetimeCells(po, pt, ne, nt, cell)
    integer, intent(in) :: po !< order of spatial elements
    integer, intent(in) :: pt !< order of temporal elements
    integer, intent(in) :: ne !< number of spatial elements
    integer, intent(in) :: nt !< number of temporal elements
    integer(C_VTK_ID), allocatable, intent(out) :: cell(:,:) !< grid cells

    integer :: c, o, i, j, l, m
    integer :: i0, j0, npe, npt

    allocate(cell(0:4, ne * nt * po * pt))

    cell(0,:) = 4 ! grid points per cell

    npe = ne*(po+1) - 1
    npt = nt*(pt+1) - 1

    c =  1 ! cell counter
    o = -1 ! cell offset

    do m = 1, nt
      j0 = (m-1)*(pt+1)
      do j = 1, pt
        do l = 1, ne
          i0 = (l-1)*(po+1)
          do i = 1, po
            cell(1,c) = o + LexicalVertexIndex(i-1 + i0, j-1 + j0, 0, npe, npt, 0)
            cell(2,c) = o + LexicalVertexIndex(i   + i0, j-1 + j0, 0, npe, npt, 0)
            cell(3,c) = o + LexicalVertexIndex(i   + i0, j   + j0, 0, npe, npt, 0)
            cell(4,c) = o + LexicalVertexIndex(i-1 + i0, j   + j0, 0, npe, npt, 0)
            c = c + 1
          end do
        end do
      end do
    end do

  end subroutine BuildLinearSpacetimeCells

  !-----------------------------------------------------------------------------
  !> Maps scalar element variables to linear grid cells

  subroutine BuildLinearScalarData(np, ns, sc, sg)
      real(RNP),      intent(in)  :: sc(:,:,:,:,:)
      real(C_DOUBLE), intent(out) :: sg(:,:)
      integer,        intent(in)  :: np, ns

      sg = reshape(sc, shape=[np, ns])

  end subroutine BuildLinearScalarData

  !=============================================================================

end module Export_VTK_Spacetime_Data__1D
