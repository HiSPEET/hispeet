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

!> summary:  Time scales of multilevel incompressible Navier-Stokes problem
!> author:   Joerg Stiller
!> date:     2025/07/28
!===============================================================================

module ML__INS__Time_Scales__3D
  use Kind_Parameters
  use XMPI
  use INS__Time_Scales__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D

  implicit none
  private

  public :: ML_INS_TimeScales_3D

  !-----------------------------------------------------------------------------
  !> Multilevel incompressible flow time scales

  type ML_INS_TimeScales_3D
    type(INS_TimeScales_3D), allocatable :: level(:)
  contains
    procedure :: Evaluate
  end type ML_INS_TimeScales_3D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of multilevel incompressible flow time scales

  subroutine Evaluate(this, ml_ins, u)
    class(ML_INS_TimeScales_3D), intent(inout) :: this
    class(ML_INS_Operator_3D), intent(in) :: ml_ins
    class(ML_MeshVariable_3D), intent(in) :: u

    real(RNP), allocatable, save :: tau_loc(:,:), tau(:,:)
    integer :: l, l_top

    ! initialization ...........................................................

    l_top = size(ml_ins % ins_op)

    if (allocated(this % level)) then
      if (size(this % level) /= l_top) then
        deallocate(this % level)
      end if
    end if

    if (.not. allocated(this % level)) then
      allocate(this % level(l_top))
    end if

    allocate(tau_loc(l_top,3))
    allocate(tau, mold = tau_loc)

    ! evaluation ...............................................................

    do l = 1, l_top
      call this % level(l) % Evaluate( ins_op = ml_ins % ins_op(l) &
                                     , u      = u % level(l) % val )
    end do

    ! globalization ............................................................

    tau_loc(:,1) = this % level % tau_conv_v
    tau_loc(:,2) = this % level % tau_conv_r
    tau_loc(:,3) = this % level % tau_diff_r

    call XMPI_Allreduce(tau_loc, tau, MPI_MIN,                  &
                        comm = ml_ins%ins_op(1)%mesh%comm_world )

    this % level % tau_conv_v = tau(:,1)
    this % level % tau_conv_r = tau(:,2)
    this % level % tau_diff_r = tau(:,3)

    ! finalization .............................................................

    deallocate(tau_loc, tau)

  end subroutine Evaluate

  !=============================================================================

end module ML__INS__Time_Scales__3D
