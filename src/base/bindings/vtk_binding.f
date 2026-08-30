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

!> summary:  Export of unstructured mesh data to VTK XML
!> author:   Joerg Stiller
!> date:     2026/08/26
!===============================================================================

module VTK_Binding
  use C_Binding
  use Execution_Control
  implicit none
  private

  public :: VTK_WriteXML_Unstructured

  ! kind parameters
  integer, parameter, public :: VTK_INT32   = C_INT32_T
  integer, parameter, public :: VTK_FLOAT64 = C_DOUBLE

  ! supported cell types
  integer(VTK_INT32), parameter, public :: VTK_POLY_LINE               =  4
  integer(VTK_INT32), parameter, public :: VTK_QUAD                    =  9
  integer(VTK_INT32), parameter, public :: VTK_HEXAHEDRON              = 12
  integer(VTK_INT32), parameter, public :: VTK_TRIQUADRATIC_HEXAHEDRON = 29

contains

  !-----------------------------------------------------------------------------
  !> Writes VTK XML unstructured mesh data with unique cell type
  !>
  !> For parallel export, `piece` must be set to the partition ID and `n_pieces`
  !> to the number of partitions. In this case, each partition writes a separate
  !> VTU file. Additionally a PVTU file is generated, which combines the pieces
  !> into a single global dataset.

  subroutine VTK_WriteXML_Unstructured( points, cells, cell_type &
                                      , pa, pa_names             &
                                      , ps, ps_names             &
                                      , pv, pv_names             &
                                      , file, piece, n_pieces    )

    real(VTK_FLOAT64),            intent(in) :: points(:,:) !< mesh points
    integer(VTK_INT32),           intent(in) :: cells(:,:)  !< mesh cells
    integer(VTK_INT32),           intent(in) :: cell_type   !< VTK cell type
    integer(VTK_INT32), optional, intent(in) :: pa(:,:)     !< attributes
    character(len=*),   optional, intent(in) :: pa_names(:) !< attribute names
    real(VTK_FLOAT64),  optional, intent(in) :: ps(:,:)     !< scalars
    character(len=*),   optional, intent(in) :: ps_names(:) !< scalar names
    real(VTK_FLOAT64),  optional, intent(in) :: pv(:,:,:)   !< vectors
    character(len=*),   optional, intent(in) :: pv_names(:) !< vector names
    character(len=*),             intent(in) :: file        !< VTK file name
    integer,            optional, intent(in) :: piece       !< piece (partition)
    integer,            optional, intent(in) :: n_pieces    !< number of pieces

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    character(len=*), parameter :: NL = new_line(' ')
    character(len=:), allocatable :: rootname
    character(len=80) :: tag

    integer :: vtk  ! output unit
    integer :: nmp  ! number of mesh points
    integer :: nmc  ! number of mesh cells
    integer :: nmd  ! number of mesh dimensions = size(points,1)
    integer :: ncp  ! number of cell points
    integer :: npa  ! number of point attributes
    integer :: nps  ! number of point scalars
    integer :: npv  ! number of point vectors

    integer(VTK_INT32) :: bytes_cc  ! bytes for cell connectivity
    integer(VTK_INT32) :: bytes_co  ! bytes for cell offsets
    integer(VTK_INT32) :: bytes_ct  ! bytes for cell types
    integer(VTK_INT32) :: bytes_p   ! bytes for points
    integer(VTK_INT32) :: bytes_pa  ! bytes for point attributes
    integer(VTK_INT32) :: bytes_ps  ! bytes for point scalars
    integer(VTK_INT32) :: bytes_pv  ! bytes for point vectors

    integer(VTK_INT32) :: offset_cc ! offset of cell connectivity
    integer(VTK_INT32) :: offset_co ! offset of cell offsets
    integer(VTK_INT32) :: offset_ct ! offset of cell types
    integer(VTK_INT32) :: offset_p  ! offset of points

    integer(VTK_INT32), allocatable :: offsets_c (:) ! offsets of cells
    integer(VTK_INT32), allocatable :: offsets_pa(:) ! offsets of point attributes
    integer(VTK_INT32), allocatable :: offsets_ps(:) ! offsets of point scalars
    integer(VTK_INT32), allocatable :: offsets_pv(:) ! offsets of point vectors

    integer :: k, o

    ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! dimensions ...............................................................

    nmd = size(points, 1)
    nmp = size(points, 2)

    ncp = size(cells, 1)
    nmc = size(cells, 2)

    if (present(pa) .and. present(pa_names)) then
      npa = size(pa_names)
    else
      npa = 0
    end if

    if (present(ps) .and. present(ps_names)) then
      nps = size(ps_names)
    else
      nps = 0
    end if

    if (present(pv) .and. present(pv_names)) then
      npv = size(pv_names)
      if (size(pv,1) /= nmd) then
        call Error( 'VTK_WriteXML_Unstructured' &
                  , 'inconsistent vector dimension' &
                  , 'VTK_Binding' )
      end if
    else
      npv = 0
    end if

    ! bytes ....................................................................

    bytes_cc = 4 * nmc * ncp      ! cell connectivity
    bytes_co = 4 * nmc            ! cell offsets
    bytes_ct = 4 * nmc            ! cell types

    bytes_p  = 8 * nmp * nmd      ! points
    bytes_pa = 4 * nmp            ! point attributes
    bytes_ps = 8 * nmp            ! point scalar
    bytes_pv = 8 * nmp * nmd      ! point vector

    ! offsets ..................................................................

    offset_p  = 0
    offset_cc = offset_p  + bytes_p  + 4
    offset_co = offset_cc + bytes_cc + 4
    offset_ct = offset_co + bytes_co + 4

    allocate(offsets_c (nmc))
    allocate(offsets_pa(npa))
    allocate(offsets_ps(nps))
    allocate(offsets_pv(npv))

    do k = 1, nmc
      offsets_c(k) = k * ncp
    end do

    o = offset_ct + bytes_ct + 4
    do k = 1, npa
      offsets_pa(k) = o + (bytes_pa + 4) * (k-1)
    end do

    o = o + (bytes_pa + 4) * npa
    do k = 1, nps
      offsets_ps(k) = o + (bytes_ps + 4) * (k-1)
    end do

    o = o + (bytes_ps + 4) * nps
    do k = 1, npv
      offsets_pv(k) = o + (bytes_ps + 4) * (k-1)
    end do

    ! tag ......................................................................

    if (present(piece)) then
      if (piece < 0) then
        return ! skip empty partition
      else
        write(tag, fmt='(A2,I0)') '_p', piece
      end if
    else
      tag = ''
    end if

    ! VTU file :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    open( newunit = vtk                           &
        , file    = trim(file)//trim(tag)//'.vtu' &
        , status  = 'replace'                     &
        , access  = 'stream'                      )

    ! head .....................................................................

    write(vtk) '<?xml version="1.0"?>' // NL

    write(vtk) '<VTKFile type="UnstructuredGrid" version="0.1" ' // &
               'byte_order="LittleEndian">' // NL

    write(vtk) '  <UnstructuredGrid>' // NL

    write(vtk) '    <Piece '                          // &
               'NumberOfPoints="' // I2C(nmp) // '" ' // &
               'NumberOfCells="'  // I2C(nmc) // '">' // NL

    ! point data specification .................................................

    write(vtk) '      <PointData>' // NL

    ! attributes
    do k = 1, npa
      write(vtk) '        '                                 // &
                 '<DataArray type="Int32" '                 // &
                 'Name="' // trim(pa_names(k)) // '" '      // &
                 'format="appended" '                       // &
                 'offset="' // I2C(offsets_pa(k)) // '" />' // NL
    end do

    ! scalars
    do k = 1, nps
      write(vtk) '        '                                 // &
                 '<DataArray type="Float64" '               // &
                 'Name="' // trim(ps_names(k)) // '" '      // &
                 'format="appended" '                       // &
                 'offset="' // I2C(offsets_ps(k)) // '" />' // NL
    end do

    ! vectors
    do k = 1, npv
      write(vtk) '        '                                  // &
                 '<DataArray type="Float64" '                // &
                 'Name="' // trim(pv_names(k)) // '" '       // &
                 'NumberOfComponents="' // I2C(nmd) // '" '  // &
                 'format="appended" '                        // &
                 'offset="' // I2C(offsets_pv(k)) // '" />'  // NL
    end do

    write(vtk) '      </PointData>' // NL

    ! mesh specification .......................................................

    ! points
    write(vtk) '      <Points>' // NL
    write(vtk) '        <DataArray type="Float64" '                  // &
                          'NumberOfComponents="' // I2C(nmd) // '" ' // &
                          'format="appended" '                       // &
                          'offset="' // I2C(offset_p) // '" />'      // NL
    write(vtk) '      </Points>' // NL

    ! cells
    write(vtk) '      <Cells>' // NL
    write(vtk) '        <DataArray type="Int32" '                // &
                          'Name="connectivity" '                 // &
                          'format="appended" '                   // &
                          'offset="' // I2C(offset_cc) // '" />' // NL
    write(vtk) '        <DataArray type="Int32" '                // &
                          'Name="offsets" '                      // &
                          'format="appended" '                   // &
                          'offset="' // I2C(offset_co) // '" />' // NL
    write(vtk) '        <DataArray type="Int32" '                // &
                          'Name="types" '                        // &
                          'format="appended" '                   // &
                          'offset="' // I2C(offset_ct) // '" />' // NL
    write(vtk) '      </Cells>' // NL

    write(vtk) '    </Piece>' // NL
    write(vtk) '  </UnstructuredGrid>' // NL

    ! appended data ............................................................

    write(vtk) '  <AppendedData encoding="raw">' // NL
    write(vtk) '_'

    write(vtk) bytes_p
    write(vtk) points

    write(vtk) bytes_cc
    write(vtk) cells

    write(vtk) bytes_co
    write(vtk) offsets_c

    write(vtk) bytes_ct
    write(vtk) spread(cell_type, dim = 1, ncopies = nmc)

    do k = 1, npa
      write(vtk) bytes_pa
      write(vtk) pa(:,k)
    end do

    do k = 1, nps
      write(vtk) bytes_ps
      write(vtk) ps(:,k)
    end do

    do k = 1, npv
      write(vtk) bytes_pv
      write(vtk) pv(:,:,k)
    end do

    write(vtk) NL // '  </AppendedData>' // NL

    ! close VTU file ...........................................................

    write(vtk) '</VTKFile>' // NL
    close(vtk)

    ! PVTU file ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(piece) .and. present(n_pieces)) then
      if (piece == 0) then

        open(newunit = vtk, file = trim(file)//'.pvtu', status = 'replace')

        ! header
        write(vtk,'(1A)') '<?xml version="1.0"?>'
        write(vtk,'(3A)') '<VTKFile type="PUnstructuredGrid" version="0.1" ', &
                           'byte_order="LittleEndian">'

        write(vtk,'(2X,A)') '<PUnstructuredGrid GhostLevel="0">'

        ! start point data section
        write(vtk,'(4X,A)') '<PPointData>'

        ! point attributes
        do k = 1, npa
          write(vtk,'(6X,4A)') '<PDataArray type="Int32" ', &
                                'Name="', trim(pa_names(k)), '"/>'
        end do

        ! scalars
        do k = 1, nps
          write(vtk,'(6X,4A)') '<PDataArray type="Float64" ', &
                                'Name="', trim(ps_names(k)), '"/>'
        end do

        ! vectors
        do k = 1, npv
          write(vtk,'(6X,5A)') '<PDataArray type="Float64" ', &
                                'Name="', trim(pv_names(k)),'" ', &
                                'NumberOfComponents="3"/>'
        end do

        ! close point data section
        write(vtk,'(4X,A)') '</PPointData>'

        ! points section
        write(vtk,'(4X,A)') '<PPoints>'
        write(vtk,'(6X,A)') '<PDataArray type="Float64" NumberOfComponents="3"/>'
        write(vtk,'(4X,A)') '</PPoints>'

        ! piece sources
        rootname = trim(file(scan(file,'/\',back=.true.)+1:))

        do k = 0, n_pieces-1
          write(vtk,'(4X,3A,I0,A)') '<Piece Source="',rootname,'_p',k,'.vtu"/>'
        end do

        ! trailer
        write(vtk,'(2X,1A)') '</PUnstructuredGrid>'
        write(vtk,'(A)') '</VTKFile>'

        ! close file
        close(vtk)

      end if
    end if

  end subroutine VTK_WriteXML_Unstructured

  !-----------------------------------------------------------------------------
  !> Converts integer to character

  function I2C(i) result(c)
    integer, intent(in) :: i
    character(len=:), allocatable :: c
    character(len=16) :: buf
    write(buf,'(I16)') i
    c = trim(adjustl(buf))
  end function I2C

  !=============================================================================

end module VTK_Binding
