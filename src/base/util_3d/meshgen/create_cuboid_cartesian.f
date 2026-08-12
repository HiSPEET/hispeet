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

!> summary:  Creation of a cuboid domain with regular Cartesian mesh
!> author:   Joerg Stiller
!> date:     2021/09/29
!===============================================================================

module Create_Cuboid_Cartesian
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generate_Regular_Mesh__3D
  implicit none
  private

  public :: CreateCuboidCartesian

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cuboid domain with regular Cartesian mesh

  subroutine CreateCuboidCartesian(comm, input, mesh)

    ! arguments ................................................................

    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    type(Mesh_3D),    intent(out) :: mesh  !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: xo(3) = 0             ! corner closest to -infinity
    real(RNP) :: lx(3) = 2*PI          ! domain extensions
    integer   :: np(3) = 1             ! number of partitions per direction
    integer   :: ep(3) = 2             ! elements per partition and direction
    integer   :: pg    = 1             ! polynomial order of geometry

    namelist/cuboid_cartesian_prm/ xo, lx, np, ep, pg

    logical   :: periodic(3) = .false. ! switch for non/periodic directions
    logical   :: regular     = .true.  ! define mesh to be (not) regular
    logical   :: structured  = .true.  ! define mesh to be (un)structured

    namelist/cuboid_cartesian_prm/ periodic, regular, structured

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io
    real(RNP) :: dx(3)         ! element spacing in directions 1:3

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_cartesian_prm)
      close(io)
    end if

    call XMPI_Bcast(xo, 0, comm)
    call XMPI_Bcast(lx, 0, comm)
    call XMPI_Bcast(np, 0, comm)
    call XMPI_Bcast(ep, 0, comm)
    call XMPI_Bcast(pg, 0, comm)

    call XMPI_Bcast(periodic  , 0, comm)
    call XMPI_Bcast(regular   , 0, comm)
    call XMPI_Bcast(structured, 0, comm)

    dx = lx / (np * ep)

    ! create mesh ..............................................................

    call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)

    mesh % regular = regular
    mesh % structured = structured

  end subroutine CreateCuboidCartesian

  !=============================================================================

end module Create_Cuboid_Cartesian
