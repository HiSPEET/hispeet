!> summary:  Collection of test cases for incompressible flow
!> author:   Joerg Stiller
!> date:     2023/03/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Problem__Test_Suite__3D

  use INS__Problem__3D                    , only: INS_Problem_3D
  use INS__Problem__Poiseuille__3D        , only: INS_Problem_Poiseuille_3D
  use INS__Problem__Stokes_DKM__3D        , only: INS_Problem_Stokes_DKM_3D
  use INS__Problem__Stokes_Linke__3D      , only: INS_Problem_Stokes_Linke_3D
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
    case('Poiseuille')
      allocate(INS_Problem_Poiseuille_3D        :: problem)
    case('Stokes_DKM')
      allocate(INS_Problem_Stokes_DKM_3D        :: problem)
    case('Stokes_Linke')
      allocate(INS_Problem_Stokes_Linke_3D      :: problem)
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
