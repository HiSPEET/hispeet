!> summary:  Element operators for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2016/03/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for IP/DG-SEM
!===============================================================================

module IP_Element_Operators_1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Standard_Operators_1D
  implicit none
  private

  public :: IP_ElementOperators1D

  !-----------------------------------------------------------------------------
  !> Element operators for symmetric interior penalty IP/DG-SEM

  type, extends(StandardOperators1D) :: IP_ElementOperators1D

    real(RNP) :: penalty = 2       !< penalty parameter > 1
    logical   :: hybrid  = .false. !< switch to IP-H

  contains

    generic :: New => New_IP_ElementOperators1D
    procedure, private :: New_IP_ElementOperators1D

    generic :: PenaltyFactor => PenaltyFactor_NE, PenaltyFactor_EQ
    procedure, private :: PenaltyFactor_NE
    procedure, private :: PenaltyFactor_EQ

    procedure :: GetStiffnessMatrix

  end type IP_ElementOperators1D

contains

!-------------------------------------------------------------------------------
!> Specific initialization, only required to override penalty

subroutine New_IP_ElementOperators1D(this, po, penalty, hybrid)
  class(IP_ElementOperators1D), intent(inout) :: this
  integer,           intent(in) :: po      !< polynomial order
  real(RNP),         intent(in) :: penalty !< penalty parameter > 1 [2]
  logical, optional, intent(in) :: hybrid  !< switch to IP-H

  ! standard operators
  call this%New(po)

  this % penalty = penalty
  if (present(hybrid)) then
    this % hybrid = hybrid
  end if

end subroutine New_IP_ElementOperators1D

!-------------------------------------------------------------------------------
!> Penalty factor for non-equidistant spacing

real(RNP) function PenaltyFactor_NE(this, dx) result(mu)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in) :: dx(2)  !< element extensions

  mu = this%penalty/4 * this%po * (this%po + 1) * (1/dx(1) + 1/dx(2))

end function PenaltyFactor_NE

!-------------------------------------------------------------------------------
!> Penalty factor for equidistant spacing

real(RNP) function PenaltyFactor_EQ(this, dx) result(mu)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in) :: dx  !< element extension

  mu = this%penalty/4 * this%po * (this%po + 1) * 2/dx

end function PenaltyFactor_EQ

!-------------------------------------------------------------------------------
!> Returns the 1D element stiffness matrix for the interior penalty DGM
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:P,0:P,-1:1)`, where `P = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine GetStiffnessMatrix(this, dx, bc, Le)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  integer   :: P, i, j
  real(RNP) :: g(-1:1), mu_0, mu_P, c_0, c_P, h_0, h_P
  real(RNP), allocatable :: delta_0(:), delta_P(:)

  ! initialization .............................................................

  P = this % po
  g = ONE / dx

  mu_0 = this % PenaltyFactor(dx(-1:0))
  mu_P = this % PenaltyFactor(dx( 0:1))

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

      if (this%hybrid) then
        h_0 = 1 / (dx(-1) * dx(0) * mu_0)
        do j = 0, P
        do i = 0, P
          Le(i,j,-1) = Le(i,j,-1) + h_0 * Ds(0,i) * Ds(P,j)
        end do
        end do
      end if

    end if

    ! own contribution (Le⁰) ...................................................

    do j = 0, P
    do i = 0, P
      Le(i,j,0) = 2 * g(0) * Ls(i,j)                           &

                + c_0 * (   g(0) * Ds   (0,i) * delta_0(j)    &
                          + g(0) * delta_0(i) * Ds   (0,j)    &
                          + mu_0 * delta_0(i) * delta_0(j) )  &

                + c_P * ( - g(0) * Ds   (P,i) * delta_P(j)    &
                          - g(0) * delta_P(i) * Ds   (P,j)    &
                          + mu_P * delta_P(i) * delta_P(j) )
    end do
    end do

    if (this%hybrid) then
      h_0 = 1 / (dx(0) * dx(0) * mu_0)
      h_P = 1 / (dx(0) * dx(0) * mu_P)
      do j = 0, P
      do i = 0, P
        Le(i,j,0) = Le(i,j,0) + h_0 * Ds(0,i) * Ds(0,j)  &
                              + h_P * Ds(P,i) * Ds(P,j)
      end do
      end do
    end if

    ! contribution from following element (Le⁺) ................................

    if (scan(bc(2), 'DN') > 0) then

      Le(:,:, 1) = 0

    else

      do j = 0, P
      do i = 0, P
        Le(i,j,1) =   g(0) * Ds   (P,i) * delta_0(j)  &
                    - g(1) * delta_P(i) * Ds   (0,j)  &
                    - mu_P * delta_P(i) * delta_0(j)
      end do
      end do

      if (this%hybrid) then
        h_P = 1 / (dx(0) * dx(1) * mu_P)
        do j = 0, P
        do i = 0, P
          Le(i,j,1) = Le(i,j,1) + h_P * Ds(P,i) * Ds(0,j)
        end do
        end do
      end if

    end if

  end associate

end subroutine GetStiffnessMatrix

!===============================================================================

end module IP_Element_Operators_1D
