!> summary:  1D Burgers breaking wave problem
!> author:   Joerg Stiller
!> date:     2019/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Burgers_Problem_1D__Breaking_Wave

  use Kind_Parameters,   only: RNP
  use Constants,         only: HALF, PI
  use Execution_Control
  use Burgers_Problem_1D

  implicit none
  private

  public :: BurgersProblem1D_BreakingWave

  !-----------------------------------------------------------------------------
  !> Type 1D Burgers breaking wave problem

  type, extends(BurgersProblem1D) :: BurgersProblem1D_BreakingWave
  contains
    procedure :: SetProblem
    procedure :: InitialValues
  end type BurgersProblem1D_BreakingWave

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the Burgers problem

  subroutine SetProblem(problem, file)
    class(BurgersProblem1D_BreakingWave), intent(inout) :: problem
    character(len=*), optional, intent(in) :: file !< input file (*.prm)

    real(RNP) :: xb1    = -1  ! position of first boundary  (left)
    real(RNP) :: xb2    =  1  ! position of second boundary (right)
    character :: bc(2)  = 'P' ! boundary conditions
    real(RNP) :: nu_0r  =  0  ! constant regular  viscosity
    real(RNP) :: nu_0s  =  0  ! constant spectral viscosity

    namelist /Burgers1D_BreakingWave_Parameters/ nu_0r, nu_0s

    logical :: exists, opened
    integer :: prm

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Warning('SetProblem',                                      &
                     'Input file "'//trim(file)//'.prm" not found, '//  &
                     'using defaults',                                  &
                     'Burgers_Problem_1D__Breaking_Wave')
      end if

    else
      exists = .false.
    end if

    ! read parameters ..........................................................

    if (exists) then
      read(prm, nml=Burgers1D_BreakingWave_Parameters)
      if (.not. opened) close(prm)
    end if

    ! set parameters ...........................................................

    problem % nc    = 1
    problem % xb1   = xb1
    problem % xb2   = xb2
    problem % bc    = reshape(bc, shape = [2,1])
    problem % nu_0r = nu_0r
    problem % nu_0s = nu_0s

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  function InitialValues(problem) result(u)
    class(BurgersProblem1D_BreakingWave), intent(in) :: problem
    real(RNP) :: u(0:problem%eop%po, problem%ne, problem%nc)

    u(:,:,1) = HALF - sin(PI * problem % x)

  end function InitialValues

  !=============================================================================

end module Burgers_Problem_1D__Breaking_Wave
