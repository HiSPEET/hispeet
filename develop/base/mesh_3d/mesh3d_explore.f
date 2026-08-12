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

!> summary:  Verification of inbuilt regular meh generation
!> author:   Joerg Stiller
!> date:     2020/12/29
!===============================================================================

program Mesh3d_Explore
  use Kind_Parameters, only: RNP
  use Constants
  use XMPI
  use Mesh_Element__3D
  use Mesh__3D
  use Generate_Regular_Mesh__3D
  use Element_Transfer_Buffer__3D
  use Verify_Mesh__3D
  use Assembly__3D
  use Export_VTK_Mesh_SFC__3D
  implicit none

  character(len=*), parameter :: input_file = 'mesh3d_explore.prm'
  real(RNP) :: xo(3) = -1                ! corner closest to -infinity
  real(RNP) :: lx(3) =  2                ! domain extensions
  integer   :: np(3) = [2,1,1]           ! number of partitions per direction
  integer   :: ep(3) = [1,1,1]           ! elements per partition and direction
  integer   :: pg    =  1                ! degree of geometry description
  integer   :: po    =  2                ! degree of mesh variables
  logical   :: periodic(3) = .false.     ! periodic directions set true
  namelist/input/ xo, lx, np, ep, pg, po, periodic

  logical   :: test_avg   = .false.      ! perform averaging test
  logical   :: export_sfc = .false.      ! export space filling curve
  namelist/control/ test_avg, export_sfc

  type(MPI_Comm) :: comm                 ! MPI communicator
  integer        :: rank                 ! local MPI rank
  integer        :: n_proc               ! number of MPI processes

  type(Mesh_3D) :: mesh         ! mesh partition
  real(RNP), allocatable :: x(:,:,:,:,:) ! mesh points
  real(RNP)              :: dx(3)        ! element spacing in directions 1:3

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), allocatable         :: v(:,:,:,:)     ! test variable
  real(RNP), pointer             :: r(:,:,:,:)     ! reference variable
  real(RNP), pointer             :: e(:,:,:,:)     ! error

  type(ElementTransferBuffer_3D), asynchronous, allocatable :: buf_v

  real(RNP) :: kappa(3), y(3), err, err_loc
  logical   :: passed, all_passed
  integer   :: io, part
  integer   :: i, j, k, l
  integer   :: i_err, j_err, k_err, l_err

  call XMPI_Init()
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  if (rank == 0) then
    open(newunit = io, file = input_file)
    read(io, nml = input)
    read(io, nml = control)
    close(io)
  end if

  call XMPI_Bcast(xo      , 0, comm)
  call XMPI_Bcast(lx      , 0, comm)
  call XMPI_Bcast(np      , 0, comm)
  call XMPI_Bcast(ep      , 0, comm)
  call XMPI_Bcast(pg      , 0, comm)
  call XMPI_Bcast(po      , 0, comm)
  call XMPI_Bcast(periodic, 0, comm)

  call XMPI_Bcast(test_avg, 0, comm)

  dx = lx / (np * ep)

  ! verification ...............................................................

  call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)
  call VerifyMesh_3D(mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)
  if (rank == 0) then
    write(*,'(/,A,G0)') 'VerifyMesh_3D: passed = ', all_passed
  end if

  ! averaging test .............................................................

  if (test_avg) then

    call mesh % GetPoints(po, 'L', x)

    allocate(v  (0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost))
    allocate(var(0:po, 0:po, 0:po, mesh%n_elem, 2))
    r(0:,0:,0:,1:) => var(:,:,:,:,1)
    e(0:,0:,0:,1:) => var(:,:,:,:,2)

    kappa = 2 * PI / lx

    do l = 1, mesh%n_elem
      do k = 0, po
      do j = 0, po
      do i = 0, po
        y(1) = x(i,j,k,l,1) - xo(1)
        y(2) = x(i,j,k,l,2) - xo(2)
        y(3) = x(i,j,k,l,3) - xo(3)
        r(i,j,k,l) = cos(sum(kappa * y))
        v(i,j,k,l) = r(i,j,k,l)
      end do
      end do
      end do
    end do

    buf_v = ElementTransferBuffer_3D(mesh, v)
    call Assembly_3D(mesh, v, buf_v, avg=.true.)

    err_loc = 0
    do l = 1, mesh%n_elem
      do k = 0, po
      do j = 0, po
      do i = 0, po
        e(i,j,k,l) = v(i,j,k,l) - r(i,j,k,l)
        if (abs(e(i,j,k,l)) > err_loc) then
          err_loc = abs(e(i,j,k,l))
          i_err   = i
          j_err   = j
          k_err   = k
          l_err   = l
        end if
      end do
      end do
      end do
    end do
    call XMPI_Reduce(err_loc, err, MPI_MAX, 0, comm)
    if (rank == 0) then
      write(*,'(/,A,ES10.3)') 'Average over element boundaries: err = ', err
    end if
    if (err_loc > epsilon(ONE)) then
      write(*,'(I5,A,ES10.3,A,4I5)') &
        rank,': err = ', err_loc,' at ', i_err, j_err, k_err, l_err
    end if
    call MPI_Barrier(comm)

  end if

  ! SFC export .................................................................

  if (export_sfc .and. mesh % has_sfc) then
    call ExportVTK_MeshSFC(mesh, file = 'mesh3d_explore_sfc')
  end if

  ! interactive exploration ....................................................

  do
    if (rank == 0) then
      write(*,'(/,A,I0,A)') 'partition: ', mesh%n_parts,' > part = '
      read(*,*) part
    end if
    call XMPI_Bcast(part, 0, comm)

    if (part < 0) exit

    if (part == rank) then
      call ShowMeshProperties(mesh)
    end if
    call MPI_Barrier(comm)

    do
      if (rank == 0) then
        write(*,'(/,A,I0,A)',advance='NO') 'element: ', mesh%n_elem,' >= l = '
        read(*,*) l
      end if
      call XMPI_Bcast(l, 0, comm)
      if (l < 1) exit
      if (part == rank) then
        if (l <= mesh % n_elem) then
          call ShowMeshElement(mesh % element(l))
        else
          call ShowMeshElement(mesh % ghost(l - mesh%n_elem))
        end if
      end if
      call MPI_Barrier(comm)
    end do

  end do

  call MPI_Finalize()

contains

  subroutine ShowMeshProperties(mesh)
    class(Mesh_3D), intent(in) :: mesh

    integer :: min_vert, max_vert
    integer :: min_edge, max_edge
    integer :: min_face, max_face
    integer :: min_elem, max_elem
    integer :: min_n_nb, max_n_nb
    integer :: i

    write(*,*)
    write(*,'(A,3(1X,G0))') 'n_vert     =', mesh % n_vert
    write(*,'(A,3(1X,G0))') 'n_edge     =', mesh % n_edge
    write(*,'(A,3(1X,G0))') 'n_face     =', mesh % n_face
    write(*,'(A,3(1X,G0))') 'n_elem     =', mesh % n_elem
    write(*,'(A,3(1X,G0))') 'n_ghost    =', mesh % n_ghost
    write(*,'(A,3(1X,G0))') 'p_geom     =', mesh % p_geom
    write(*,'(A,3(1X,G0))') 'n_bound    =', mesh % n_bound
    write(*,'(A,3(1X,G0))') 'n_part     =', mesh % n_parts
    write(*,'(A,3(1X,G0))') 'part       =', mesh % part
    write(*,'(A,3(1X,G0))') 'structured =', mesh % structured
    write(*,'(A,3(1X,G0))') 'regular    =', mesh % regular
    write(*,'(A,3(1X,G0))') 'n_elem_1   =', mesh % n_elem_1
    write(*,'(A,3(1X,G0))') 'n_elem_2   =', mesh % n_elem_2
    write(*,'(A,3(1X,G0))') 'n_elem_3   =', mesh % n_elem_3
    write(*,'(A,3(1X,G0))') 'dx         =', mesh % dx

    min_vert = huge(1);  max_vert = -huge(1)
    min_face = huge(1);  max_face = -huge(1)
    min_edge = huge(1);  max_edge = -huge(1)
    min_elem = huge(1);  max_elem = -huge(1)
    min_n_nb = huge(1);  max_n_nb = -huge(1)

    do i = 1, mesh % n_elem
      associate(element => mesh % element(i))
        min_vert = min( min_vert, minval(element % vertex % id) )
        max_vert = max( max_vert, maxval(element % vertex % id) )
        min_edge = min( min_edge, minval(element % edge   % id) )
        max_edge = max( max_edge, maxval(element % edge   % id) )
        min_face = min( min_face, minval(element % face   % id) )
        max_face = max( max_face, maxval(element % face   % id) )
        min_elem = min( min_elem, element % id )
        max_elem = max( max_elem, element % id )
        min_n_nb = min( min_n_nb, size(element % neighbor) )
        max_n_nb = max( max_n_nb, size(element % neighbor) )
      end associate
    end do
    write(*,*)
    write(*,'(A,3(1X,G0))') 'min/max_vert =', min_vert, max_vert
    write(*,'(A,3(1X,G0))') 'min/max_edge =', min_edge, max_edge
    write(*,'(A,3(1X,G0))') 'min/max_face =', min_face, max_face
    write(*,'(A,3(1X,G0))') 'min/max_elem =', min_elem, max_elem
    write(*,'(A,3(1X,G0))') 'min/max_n_nb =', min_n_nb, max_n_nb

  end subroutine ShowMeshProperties

  subroutine ShowMeshElement(element)
    class(MeshElement_3D), intent(in) :: element

    integer :: i, j, j1, j2, k

    write(*,*)
    write(*,'(A,99(1X,I5))') 'id            =', element%id
    write(*,'(A,99(1X,I5))') 'vertex % id   =', element%vertex%id
    write(*,'(A,99(1X,I5))') 'vertex % rank =', element%vertex%rank
    write(*,'(A,99(1X,I5))') 'vertex % val  =', element%vertex%val
    write(*,'(A,99(1X,I5))') 'edge   % id   =', element%edge%id
    write(*,'(A,99(1X,I5))') 'edge   % rank =', element%edge%rank
    write(*,'(A,99(1X,I5))') 'edge   % val  =', element%edge%val
    write(*,'(A,99(1X,I5))') 'face   % id   =', element%face%id
    write(*,'(A,99(1X,I5))') 'face   % rank =', element%face%rank
    write(*,'(A,99(1X,I5))') 'face   % val  =', element%face%val
    write(*,*)

    if (.not.allocated(element%neighbor)) return

    write(*,'(A,99(1X,I5))') 'n_neighbor  =', size(element%neighbor)
    write(*,*)
    k = 0
    do i = 1, 6
      j1 = element % face(i) % i_neighbor
      j2 = element % face(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,4I5)') k, '  face  ', i, &
            '  id,part,component,orientation =', element % neighbor(j)
      end do
    end do
    do i = 1, 12
      j1 = element % edge(i) % i_neighbor
      j2 = element % edge(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,4I5)') k, '  edge  ', i, &
            '  id,part,component,orientation =', element % neighbor(j)
      end do
    end do
    do i = 1, 8
      j1 = element % vertex(i) % i_neighbor
      j2 = element % vertex(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,4I5)') k, '  vertex', i, &
            '  id,part,component,orientation =', element % neighbor(j)
      end do
    end do

  end subroutine ShowMeshElement

end program Mesh3d_Explore
