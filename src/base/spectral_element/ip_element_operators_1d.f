!> summary:  Element operators for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2016/03/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for IP/DG-SEM
!===============================================================================

module IP_Element_Operators_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  implicit none
  private

  public :: IP_ElementOperators1D

  !-----------------------------------------------------------------------------
  !> Element operators for the symmetric interior penalty IP/DG-SEM
  !>
  !> The operator accommodates classical IP as well as hybridizable IP-H. The
  !> particular method is selected during initialization and controlled by then
  !> component `hybrid`.
  !>
  !> In case of IP-H, the operator provides  the column matrix of generalized
  !> eigenvectors `S` and the diagonal matrix of eigenvalues `Λ = Lambda` such
  !> that
  !>
  !>     Sᵀ Lᵢᵢ S = Λ
  !>     Sᵀ Mᵢᵢ S = I
  !>
  !> where `Lᵢᵢ` and `Mᵢᵢ` the standard stiffness matrix and the standard
  !> diagonal mass matrix restricted to the interior points.

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
    procedure :: GetEllipticEigensystem
    procedure :: GetEllipticSuboperators

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
!>
!> The stiffness matrix is available in two forms
!>
!>   * `'primal'`: all numeric fluxes û are eliminated (default)
!>   * `'flux'`  : û is retained, all corresponding terms are removed from `Le`
!>
!> Except for the single element case, i.e. `all(bc /= '')`, the flux form can
!> be activated by passing `form = 'flux'`.

subroutine GetStiffnessMatrix(this, dx, bc, Le, form)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N','P'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix
  character(len=*), optional, intent(in) :: form !< operator form ['primal']

  logical   :: primal
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

  if (present(form)) then
    primal = form == 'primal' .or. all(bc /= ' ')
  else
    primal = .true.
  end if

  if (primal) then

    ! left boundary condition
    select case(bc(1))
    case('D')    ! Dirichlet
      c_0 = 2
    case('N')    ! Neumann
      c_0 = 0
    case default ! none
      c_0 = 1
    end select

    ! right boundary condition
    select case(bc(2))
    case('D')    ! Dirichlet
      c_P = 2
    case('N')    ! Neumann
      c_P = 0
    case default ! none
      c_P = 1
    end select

  else

    ! left boundary condition
    select case(bc(1))
    case('N')    ! Neumann
      c_0 = 0
    case default
      c_0 = 2
    end select

    ! right boundary condition
    select case(bc(2))
    case('N')    ! Neumann
      c_P = 0
    case default
      c_P = 2
    end select

  end if

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

      if (primal .and. this%hybrid) then
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
      Le(i,j,0) = 2 * g(0) * Ls(i,j)                          &

                + c_0 * (   g(0) * Ds   (0,i) * delta_0(j)    &
                          + g(0) * delta_0(i) * Ds   (0,j)    &
                          + mu_0 * delta_0(i) * delta_0(j) )  &

                + c_P * ( - g(0) * Ds   (P,i) * delta_P(j)    &
                          - g(0) * delta_P(i) * Ds   (P,j)    &
                          + mu_P * delta_P(i) * delta_P(j) )
    end do
    end do

    if (primal .and. this%hybrid) then

      select case(bc(1))
      case(' ','P')
        h_0 = 1 / (dx(0) * dx(0) * mu_0)
        do j = 0, P
        do i = 0, P
          Le(i,j,0) = Le(i,j,0) - h_0 * Ds(0,i) * Ds(0,j)
        end do
        end do
      end select

      select case(bc(2))
      case(' ','P')
        h_P = 1 / (dx(0) * dx(0) * mu_P)
        do j = 0, P
        do i = 0, P
          Le(i,j,0) = Le(i,j,0) - h_P * Ds(P,i) * Ds(P,j)
        end do
        end do
      end select

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

      if (primal .and. this%hybrid) then
        h_P = 1 / (dx(0) * dx(1) * mu_P)
        do j = 0, P
        do i = 0, P
          Le(i,j,1) = Le(i,j,1) + h_P * Ds(P,i) * Ds(0,j)
        end do
        end do
      end if

    end if

    ! special case: single periodic element ....................................

    if (all(bc == 'P')) then
      Le(:,:, 0) = Le(:,:,0) + Le(:,:,-1) + Le(:,:,1)
      Le(:,:,-1) = 0
      Le(:,:, 1) = 0
    end if

  end associate

end subroutine GetStiffnessMatrix

!-------------------------------------------------------------------------------
!> Provides the generalized eigensystem for interior stiffness and mass matrices
!>
!> Returns the column matrix of generalized eigenvectors `S` and the diagonal
!> matrix of eigenvalues `Λ = Lambda` to the interior element stiffness matrix
!> `Lᵢᵢ` and diagonal mass matrix `Mᵢᵢ` of the hybridized element system such
!> that
!>
!>     Sᵀ Lᵢᵢ S = Λ
!>     Sᵀ Mᵢᵢ S = I

subroutine GetEllipticEigensystem(this, dx, bc, S, Lambda)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)   !< element extensions
  character, intent(in)  :: bc(2)      !< boundary conditions {'','D','N','P'}
  real(RNP), intent(out) :: S(0:,0:)   !< eigenvectors
  real(RNP), intent(out) :: Lambda(0:) !< eigenvalues

  real(RNP), parameter   :: eps = epsilon(1.0)
  real(RNP), allocatable :: Mii(:), Lii(:,:,:)

  if (.not. this%hybrid) then
    call Error( 'GetEllipticEigensystem', &
                'available only for hybridizable IP', &
                'IP_Element_Operators_1D' )
  end if

  ! hybrid element operators
  allocate(Mii(0:this%po), Lii(0:this%po, 0:this%po, -1:1))
  Mii = dx(0)/2 * this % w
  call this % GetStiffnessMatrix(dx, bc, Lii, form='flux')

  ! solve eigenproblem
  call SolveGeneralizedEigenproblem(Lii(:,:,0), Mii, Lambda, S)

  ! singular case
  if (all(bc == 'N') .or. all(bc == 'P')) then
    Lambda(0) = 0
  end if

end subroutine GetEllipticEigensystem

!-------------------------------------------------------------------------------
!> Computes operators for hybrid IP/DG-SEM diffusion problem
!>
!> To cope with the singular case (Neumann or periodic with c = 0),
!> we use the Moore-Penrose inverse, i.e. `Aii_inv = Ã⁺`

subroutine GetEllipticSuboperators(this, dx, bc, c, nu, Aib, Aii_inv)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)       !< element extensions
  character, intent(in)  :: bc(2)          !< boundary conds {'','D','N','P'}
  real(RNP), intent(in)  :: c              !< coefficient of linear term
  real(RNP), intent(in)  :: nu             !< diffusivity
  real(RNP), intent(out) :: Aib(0:,:)      !< interior-boundary part, Â(0:P,1:2)
  real(RNP), intent(out) :: Aii_inv(0:,0:) !< inv interior part, Ã⁺(0:P,0:P)

  real(RNP), allocatable :: S(:,:), Lambda(:), D_inv(:)
  real(RNP), allocatable :: delta_0(:), delta_P(:)
  real(RNP) :: mu_0, mu_P
  integer   :: P, i, j

  ! initialization .............................................................

  P = this % po

  mu_0 = this % PenaltyFactor(dx(-1:0))
  mu_P = this % PenaltyFactor(dx( 0:1))

  allocate(delta_0(0:P), source = ZERO)
  delta_0(0) = ONE

  allocate(delta_P(0:P), source = ZERO)
  delta_P(P) = ONE

  allocate(S(0:P,0:P), Lambda(0:P), D_inv(0:P))
  call this % GetEllipticEigensystem(dx, bc, S, Lambda)

  ! interior-boundary part .....................................................

  associate(Ds => this%D)

    if (all(bc == 'P')) then
      Aib(:,1) =  0
      Aib(:,2) =  0
    else

      select case (bc(1))
      case('D','N')
        Aib(:,1) =  0
      case default ! interior od periodic
        Aib(:,1) = -2/dx(0) * nu * Ds(0,:)  -  2 * nu * mu_0 * delta_0
      end select

      select case (bc(2))
      case('D','N')
        Aib(:,2) =  0
      case default ! interior od periodic
        Aib(:,2) =  2/dx(0) * nu * Ds(P,:)  -  2 * nu * mu_P * delta_P
      end select

    end if

  end associate

  ! inverse interior part ......................................................

  D_inv = c + nu * Lambda
  where(abs(D_inv) > epsilon(ONE))
    D_inv = 1 / D_inv
  elsewhere
    D_inv = ZERO
  end where

  do j = 0, P
  do i = 0, P
    Aii_inv(i,j) = sum(S(i,:) * D_inv * S(j,:))
  end do
  end do

end subroutine GetEllipticSuboperators

!===============================================================================

end module IP_Element_Operators_1D
