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

!> summary:  Collection of test cases for incompressible flow
!> author:   Joerg Stiller
!> date:     2023/03/12
!===============================================================================

module INS__Problem__Test_Suite__3D

  use INS__Problem__3D                    , only: INS_Problem_3D
  use INS__Problem__Channel__3D           , only: INS_Problem_Channel_3D
  use INS__Problem__Cylinder_2D__3D       , only: INS_Problem_Cylinder2D_3D
  use INS__Problem__Hagen_Poiseuille__3D  , only: INS_Problem_HagenPoiseuille_3D
  use INS__Problem__Linear_Cascade__3D    , only: INS_Problem_LinearCascade_3D
  use INS__Problem__No_Flow__3D           , only: INS_Problem_NoFlow_3D
  use INS__Problem__Poiseuille__3D        , only: INS_Problem_Poiseuille_3D
  use INS__Problem__Stokes_DKM__3D        , only: INS_Problem_Stokes_DKM_3D
  use INS__Problem__Stokes_Linke__3D      , only: INS_Problem_Stokes_Linke_3D
  use INS__Problem__Transition_TG__3D     , only: INS_Problem_TransitionTG_3D
  use INS__Problem__Variable_Viscosity__3D, only: INS_Problem_VariableViscosity_3D
  use INS__Problem__Vortex_TG__3D         , only: INS_Problem_Vortex_TG_3D

  use XMPI

  implicit none
  private

  public :: INS_Problem_3D
  public :: Set_INS_TestProblem_3D

contains

  !-----------------------------------------------------------------------------
  !>

  subroutine Set_INS_TestProblem_3D(problem, name, file, n_bound, comm)
    class(INS_Problem_3D), allocatable, intent(out) :: problem
    character(len=*), intent(in) :: name    !< predefined problem name
    character(len=*), intent(in) :: file    !< problem input file
    integer,          intent(in) :: n_bound !< number of boundaries
    type(MPI_Comm),   intent(in) :: comm    !< MPI world communicator

    select case(name)
    case('Channel')
      allocate(INS_Problem_Channel_3D           :: problem)
    case('Cylinder2D')
      allocate(INS_Problem_Cylinder2D_3D        :: problem)
    case('HagenPoiseuille')
      allocate(INS_Problem_HagenPoiseuille_3D   :: problem)
    case('LinearCascade')
      allocate(INS_Problem_LinearCascade_3D     :: problem)
    case('NoFlow')
      allocate(INS_Problem_NoFlow_3D            :: problem)
    case('Poiseuille')
      allocate(INS_Problem_Poiseuille_3D        :: problem)
    case('Stokes_DKM')
      allocate(INS_Problem_Stokes_DKM_3D        :: problem)
    case('Stokes_Linke')
      allocate(INS_Problem_Stokes_Linke_3D      :: problem)
    case('Transition_TG')
      allocate(INS_Problem_TransitionTG_3D      :: problem)
    case('VariableViscosity')
      allocate(INS_Problem_VariableViscosity_3D :: problem)
    case('Vortex_TG')
      allocate(INS_Problem_Vortex_TG_3D         :: problem)
    case default
    end select

    call problem % SetProblem(n_bound, file, comm)

  end subroutine Set_INS_TestProblem_3D

  !=============================================================================

end module INS__Problem__Test_Suite__3D
