!> summary:  Polynomial multigrid problem setup
!> author:   Joerg Stiller
!> date:     2019/02/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid problem setup
!===============================================================================

submodule(CART__Elliptic_PMG) MP_SetProblem
  implicit none

contains

!-------------------------------------------------------------------------------
!> Initialize problem with constant isotropic diffusivity

module subroutine SetProblem_CI(this, lambda, nu, bc)
  class(PMG_Method3D), intent(inout) :: this
  real(RNP),           intent(in)    :: lambda  !< Helmholtz parameter
  real(RNP),           intent(in)    :: nu      !< diffusivity
  character,           intent(in)    :: bc(:)   !< boundary conditions

  integer :: l, l_top

  associate(level => this % level)

    l_top = ubound(level,1)

    call level(l_top) % SetTopLevelProblem_CI(lambda, nu, bc)
    do l = l_top-1, 0, -1
      call level(l) % SetCoarseProblem(fine = level(l+1))
    end do

  end associate

end subroutine SetProblem_CI

!-------------------------------------------------------------------------------
!> Initialize problem with a combination of a constant isotropic diffusivity and
!> a constant isotropic spectral diffusivity

module subroutine SetProblem_CI_svv(this, lambda, nu, nu_svv, bc)
  class(PMG_Method3D), intent(inout) :: this
  real(RNP),           intent(in)    :: lambda  !< Helmholtz parameter
  real(RNP),           intent(in)    :: nu      !< diffusivity
  real(RNP),           intent(in)    :: nu_svv  !< spectral diffusivity
  character,           intent(in)    :: bc(:)   !< boundary conditions

  integer :: l, l_top

  associate(level => this % level)

    l_top = ubound(level,1)

    call level(l_top) % SetTopLevelProblem_CI_svv(lambda, nu, nu_svv, bc)
    do l = l_top-1, 0, -1
      call level(l) % SetCoarseProblem(fine = level(l+1))
    end do

  end associate

end subroutine SetProblem_CI_svv

!-------------------------------------------------------------------------------
!> Initialize problem with variable isotropic diffusivity

module subroutine SetProblem_VI(this, lambda, nu, bc)
  class(PMG_Method3D), intent(inout) :: this
  real(RNP),           intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP),           intent(in)    :: nu(0:,0:,0:,:) !< diffusivity
  character,           intent(in)    :: bc(:)          !< boundary conditions

  integer :: l, l_top

  associate(level => this % level)

    l_top = ubound(level,1)

    call level(l_top) % SetTopLevelProblem_VI(lambda, nu, bc)
    do l = l_top-1, 0, -1
      call level(l) % SetCoarseProblem(fine = level(l+1))
    end do

  end associate

end subroutine SetProblem_VI

!===============================================================================

end submodule MP_SetProblem
