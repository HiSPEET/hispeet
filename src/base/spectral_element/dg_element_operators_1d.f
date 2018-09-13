!> summary:  Element operators for IP-DGM
!> author:   Joerg Stiller
!> date:     2016/03/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for IP-DGM
!===============================================================================

module DG_Element_Operators_1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Standard_Operators_1D
  implicit none
  private

!-------------------------------------------------------------------------------
!> Element operators for symmetric interior penalty DGM

  type, extends(StandardOperators1D), public :: DG_ElementOperators1D
  contains
    procedure :: PenaltyFactor
    procedure :: GetStiffnessMatrix
  end type DG_ElementOperators1D

contains

!-------------------------------------------------------------------------------
!> Penalty factor

real(RNP) function PenaltyFactor(this, dx, penalty) result(mu)
  class(DG_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in) :: dx(2)    !< element extensions
  real(RNP), intent(in) :: penalty  !< penalty parameter \( \mu_\star \)

  mu = penalty/4 * this%po * (this%po + 1) * (1/dx(1) + 1/dx(2))

end function PenaltyFactor

!-------------------------------------------------------------------------------
!> Returns the 1D element stiffness matrix
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:po,0:po,-1:1)`, where `po = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine GetStiffnessMatrix(this, dx, penalty, bc, Le)
  class(DG_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
  real(RNP), intent(in)  :: penalty       !< penalty parameter (> 1)
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  integer   :: P, i, j
  real(RNP) :: g(-1:1), mu_0, mu_P, c_0, c_P
  real(RNP), allocatable :: delta_0(:), delta_P(:)

  ! initialization .............................................................

  P = this % PolynomialOrder()
  g = ONE / dx

  mu_0 = this % PenaltyFactor(dx(-1:0), penalty)
  mu_P = this % PenaltyFactor(dx( 0:1), penalty)

  allocate(delta_0(0:P), source = ZERO)
  delta_0(0) = ONE

  allocate(delta_P(0:P), source = ZERO)
  delta_P(P) = ONE

  ! left boundary condition
  select case(bc(1))
  case('D')    ! Dirichlet
    c_0 = 2
  case('N')    ! Neumann
    c_0 = 0
  case default ! none
    c_0 = 1
  end select

  ! right boundary
  select case(bc(2))
  case('D')    ! Dirichlet
    c_P = 2
  case('N')    ! Neumann
    c_P = 0
  case default ! none
    c_P = 1
  end select

  associate( Ms => this%w, Ds => this%D, Ls => this%L )

    ! contribution from preceding element (Le⁻) ................................

    if (scan(bc(1), 'DN') > 0) then
      Le(:,:,-1) = 0
    else
      do j = 0, P
      do i = 0, P
        Le(i,j,-1) = - g( 0) * Ds   (0,i) * delta_P(j)  &
                     + g(-1) * delta_0(i) * Ds   (P,j)  &
                     - mu_0  * delta_0(i) * delta_P(j)
      end do
      end do
    end if

    ! own contribution (Le⁰) ...................................................

    do j = 0, P
    do i = 0, P
      Le(i,j, 0) = 2 * g(0) * Ls(i,j)                           &

                 + c_0 * (   g( 0) * Ds   (0,i) * delta_0(j)    &
                           + g( 0) * delta_0(i) * Ds   (0,j)    &
                           + mu_0  * delta_0(i) * delta_0(j) )  &

                 + c_P * ( - g( 0) * Ds   (P,i) * delta_P(j)    &
                           - g( 0) * delta_P(i) * Ds   (P,j)    &
                           + mu_P  * delta_P(i) * delta_P(j) )
    end do
    end do

    ! contribution from following element (Le⁺) ................................

    if (scan(bc(2), 'DN') > 0) then
      Le(:,:, 1) = 0
    else
      do j = 0, P
      do i = 0, P
        Le(i,j, 1) =   g( 0) * Ds   (P,i) * delta_0(j)  &
                     - g( 1) * delta_P(i) * Ds   (0,j)  &
                     - mu_P  * delta_P(i) * delta_0(j)
      end do
      end do
    end if

  end associate

end subroutine GetStiffnessMatrix

!===============================================================================

end module DG_Element_Operators_1D
