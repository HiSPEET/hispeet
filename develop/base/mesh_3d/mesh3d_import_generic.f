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

!> summary:  Test of generic mesh generation and import
!> author:   Joerg Stiller
!> date:     2021/01/09
!===============================================================================

program Mesh3d_Import_Generic
  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Element_Transfer_Buffer__3D
  use Verify_Mesh__3D
  use Assembly__3D
  use Export_VTK_Volume_Data__3D
  use Export_VTK_Mesh_SFC__3D
  implicit none

  character(len=*), parameter :: input_file = 'mesh3d_import_generic.prm'

  integer :: config = 1 ! configuration (1 cylinder, 2 annular gap, 3 diamonds)
  logical :: test_sem    = .false. ! test SEM finctionality
  logical :: test_avg    = .false. ! perform averaging test
  logical :: export_vtk  = .false. ! generate VTK file
  namelist/control/ config, test_sem, test_avg, export_vtk

  real(RNP) :: r0   = 0.5 ! inner radius  (annular gap only)
  real(RNP) :: r1   = 1   ! outer radius
  real(RNP) :: h    = 2   ! height = axial extension
  real(RNP) :: lx   = 1   ! domain length in x-direction
  real(RNP) :: ly   = 1   ! domain length in y-direction
  real(RNP) :: lz   = 1   ! domain length in z-direction
  integer   :: nx   = 2   ! num cells in x-direction
  integer   :: ny   = 2   ! num cells in y-direction
  integer   :: nz   = 2   ! num elements in z-/axial  direction
  integer   :: nr   = 2   ! num elements in radial    direction
  integer   :: np   = 4   ! num elements in azimuthal direction ≥ 3 (gap only)
  integer   :: po   = 3   ! polynomial order of mesh elements
  logical   :: periodic = .false. ! F/T for non/periodic z-/axial direction

  namelist/cylinder/ nr, nz, po, periodic
  namelist/annulus/  r0, r1, h, nr, np, nz, po, periodic
  namelist/diamonds/ lx, ly, lz, nx, ny, nz, po, periodic

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh_3D) :: generic_mesh
  type(Mesh_3D) :: mesh
  type(SpectralElementMesh_3D) :: se_mesh

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), allocatable         :: x(:,:,:,:,:)
  real(RNP), allocatable         :: v(:,:,:,:)     ! test variable
  real(RNP), pointer             :: r(:,:,:,:)     ! reference variable
  real(RNP), pointer             :: e(:,:,:,:)     ! error

  type(ElementTransferBuffer_3D), asynchronous, allocatable :: buf_v

  real(RNP), allocatable :: area(:)
  real(RNP) :: vol
  real(RNP) :: kappa(3), y(3), err
  logical   :: passed = .false.
  integer   :: io
  integer   :: i, j, k, l
  integer   :: i_err, j_err, k_err, l_err

  type(Mesh_3D) :: new_mesh
  integer :: n_proc

  call Init_MPI_Binding()
  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then

    ! read parameters ..........................................................

    open(newunit = io, file = input_file)
    read(io, nml = control)
    select case(config)
    case(1)
      read(io, nml = cylinder)
    case(2)
      read(io, nml = annulus)
    case(3)
      read(io, nml = diamonds)
    end select
    close(io)

    ! create generic mesh ......................................................

    select case(config)
    case(1)
      call generic_mesh % CreateCylinder(r1, h, nr, nz, po, periodic)
    case(2)
      call generic_mesh % CreateAnnulus(r0, r1, h, nr, np, nz, po, periodic)
    case(3)
      call generic_mesh % CreateDiamonds(lx, ly, lz, nx, ny, nz, po, periodic)
    case default
    end select

    ! import generic mesh ......................................................

    call mesh % ImportGenericMesh(generic_mesh, comm = comm)
    call VerifyMesh_3D(mesh, passed)
    write(*,'(/,A,G0,/)') 'verification of imported mesh: passed = ', passed

    ! spectral element functionality ...........................................

    if (test_sem) then

      se_mesh = SpectralElementMesh_3D(mesh, po)

      allocate(area(mesh % n_bound))
      call se_mesh % Get_Volume(vol)
      call se_mesh % Get_SurfaceAreas(area)

      write(*,'(A)') 'Spectral element mesh'
      write(*,'(2X,A,G0)') 'volume  = ', vol
      do i = 1, mesh % n_bound
        write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
      end do

    end if

    ! set up data ..............................................................

    if (test_avg .or. export_vtk) then
      call mesh % GetPoints(po, 'L', x)
      allocate(v  (0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost))
      allocate(var(0:po, 0:po, 0:po, mesh%n_elem, 2), source = ZERO)
      r(0:,0:,0:,1:) => var(:,:,:,:,1)
      e(0:,0:,0:,1:) => var(:,:,:,:,2)
    end if

    ! averaging test ...........................................................

    if (test_avg) then

      kappa = 2 * PI / h

      do l = 1, mesh%n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          y(1) = x(i,j,k,l,1)
          y(2) = x(i,j,k,l,2)
          y(3) = x(i,j,k,l,3)
          r(i,j,k,l) = cos(sum(kappa * y))
          v(i,j,k,l) = r(i,j,k,l)
        end do
        end do
        end do
      end do

      buf_v = ElementTransferBuffer_3D(mesh, v)
      call Assembly_3D(mesh, v, buf_v, avg=.true.)

      err = 0
      do l = 1, mesh%n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          e(i,j,k,l) = v(i,j,k,l) - r(i,j,k,l)
          if (abs(e(i,j,k,l)) > err) then
            err   = abs(e(i,j,k,l))
            i_err = i
            j_err = j
            k_err = k
            l_err = l
          end if
        end do
        end do
        end do
      end do
      write(*,'(A,ES10.3,A,4I5)') 'Average over element boundaries: err = ', &
                                  err, ' at ', i_err, j_err, k_err, l_err
    end if

    ! export mesh and data .....................................................

    if (export_vtk) then
      call ExportVTK_VolumeData( x, var                   &
                               , sname   = ['r','e']      &
                               , file    = 'element_mesh' &
                               , part    = mesh % part    &
                               , n_parts = mesh % n_parts )

      call mesh % GetCuboids(x)
      call ExportVTK_VolumeData( x                        &
                               , file    = 'cuboid_mesh'  &
                               , part    = mesh % part    &
                               , n_parts = mesh % n_parts )

     if (mesh % has_sfc) then
       call ExportVTK_MeshSFC(mesh, file = 'sfc')
     end if

    end if

  end if

  call MPI_Finalize()

end program Mesh3d_Import_Generic
