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

!> summary:  Creation of a 3x3x3 cuboid with a rotated center element
!> author:   Joerg Stiller
!> date:     2022/03/20
!===============================================================================

module Create_Cuboid_OneRotated
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generic_Mesh__3D
  implicit none
  private

  public :: CreateCuboidOneRotated

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cylindrical domain with an unstructured mesh

  subroutine CreateCuboidOneRotated(comm, input, mesh)

    ! arguments ................................................................

    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    type(Mesh_3D),    intent(out) :: mesh  !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: xo(3) = 0             ! corner closest to -infinity
    real(RNP) :: lx(3) = 2*PI          ! domain extensions
    integer   :: pg = 3                ! polynomial order of geometry
    integer   :: rotation(3) = 0       ! x/y/z-rotation of center element
    logical   :: periodic(3) = .false. ! F/T for non/periodic directions

    namelist/cuboid_onerotated_prm/ xo, lx, pg, rotation, periodic

    ! auxiliary variables ......................................................

    integer :: rank
    integer :: io

    type(GenericMesh_3D) :: generic_mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_onerotated_prm)
      close(io)
    end if

    call XMPI_Bcast(xo      , 0, comm)
    call XMPI_Bcast(lx      , 0, comm)
    call XMPI_Bcast(pg      , 0, comm)
    call XMPI_Bcast(rotation, 0, comm)
    call XMPI_Bcast(periodic, 0, comm)

    ! create mesh ..............................................................

    call generic_mesh % CreateOneRotated(xo, lx, pg, rotation, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  end subroutine CreateCuboidOneRotated

  !=============================================================================

end module Create_Cuboid_OneRotated
