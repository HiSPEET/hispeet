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

!> summary:  Assembly of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2021/03/03
!===============================================================================

module Assembly__3D
  use Kind_Parameters, only: RNP
  use Constants      , only: ONE, ZERO
  use Execution_Control
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  implicit none
  private

  public :: Assembly_3D

  interface Assembly_3D
    module procedure :: Assembly_3D_0
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Unstructured assembly and, optionally, averaging of a scalar variable

    module subroutine Assembly_3D_0U(mesh, u, buf_u, avg)
      !> mesh partition
      class(Mesh_3D), intent(in) :: mesh
      !> mesh variable, including ghost entries
      real(RNP), intent(inout) :: u(:,:,:,:)
      !> MPI transfer buffer
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_u
      !> switch for averaging over element boundaries [F]
      logical, optional, intent(in) :: avg
    end subroutine Assembly_3D_0U

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Assembly and, optionally, averaging of a scalar variable

  subroutine Assembly_3D_0(mesh, u, buf_u, avg)
    !> mesh partition
    class(Mesh_3D), intent(in) :: mesh
    !> mesh variable, including ghost entries
    real(RNP), intent(inout) :: u(:,:,:,:)
    !> MPI transfer buffer
    type(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_u
    !> switch for averaging over element boundaries [F]
    logical, optional, intent(in) :: avg

    !$omp master
    if (size(u,4) /= mesh%n_elem + mesh%n_ghost) then
      call Error('Assembly_3D_0', 'Dimension 4 of u does not match')
    end if
    !$omp end master

    ! insert special treatment of structured case here
    call Assembly_3D_0U(mesh, u, buf_u, avg)

  end subroutine Assembly_3D_0

  !=============================================================================

end module Assembly__3D
