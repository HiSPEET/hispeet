!> summary:  Element operators for continuous Galerkin-SEM
!> author:   Joerg Stiller
!> date:     2018/09/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Element operators for continuous Galerkin-SEM
!===============================================================================

module CG_Element_Operators_1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Eigenproblems,   only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Element operators for continuous Galerkin-SEM

  type, extends(StandardOperators1D), public :: CG_ElementOperators1D
  contains
    procedure :: GetStiffnessMatrix
    procedure :: GetEllipticEigensystem
    procedure :: GetEllipticSuboperators
  end type CG_ElementOperators1D

contains

!-------------------------------------------------------------------------------
!> Returns the 1D element stiffness matrix for the continuous Galerkin SEM
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:P,0:P,-1:1)`, where `P = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine GetStiffnessMatrix(this, dx, bc, Le)
  class(CG_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  integer   :: P, i, j
  real(RNP) :: g(-1:1)
  real(RNP), allocatable :: delta_0(:), delta_P(:)

  ! initialization .............................................................

  P = this % po
  g = 2 / dx

  allocate(delta_0(0:P), source = ZERO)
  delta_0(0) = ONE

  allocate(delta_P(0:P), source = ZERO)
  delta_P(P) = ONE

  associate(Ls => this%L)

    ! contribution from preceding element (Le⁻) ................................

    if (scan(bc(1), 'DN') > 0) then
      Le(:,:,-1) = 0
    else
      do j = 0, P
      do i = 0, P
        Le(i,j,-1) = g(-1) * delta_0(i) * Ls(P,j)
      end do
      end do
    end if

    ! own contribution (Le⁰) ...................................................

    Le(:,:,0) = g(0) * Ls

    ! nullify Dirichlet entries
    if (bc(1) == 'D') then
      Le(0,:,0) = 0
      Le(:,0,0) = 0
    end if
    if (bc(2) == 'D') then
      Le(P,:,0) = 0
      Le(:,P,0) = 0
    end if

    ! contribution from following element (Le⁺) ................................

    if (scan(bc(2), 'DN') > 0) then
      Le(:,:,1) = 0
    else
      do j = 0, P
      do i = 0, P
        Le(i,j,1) = g(1) * delta_P(i) * Ls(0,j)
      end do
      end do
    end if

  end associate

end subroutine GetStiffnessMatrix

!-------------------------------------------------------------------------------
!> Provides the generalized eigensystem for interior stiffness and mass matrices
!>
!> Returns the column matrix of generalized eigenvectors `S` and the diagonal
!> matrix of eigenvalues `Λ = Lambda` to the interior element stiffness matrix
!> `Lᵢᵢ` and diagonal mass matrix `Mᵢᵢ` such that
!>
!>     Sᵀ Lᵢᵢ S = Λ
!>     Sᵀ Mᵢᵢ S = I

subroutine GetEllipticEigensystem(this, dx, S, Lambda)
  class(CG_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx         !< element length
  real(RNP), intent(out) :: S(:,:)     !< eigenvectors
  real(RNP), intent(out) :: Lambda(:)  !< eigenvalues

  real(RNP), allocatable :: Mii(:), Lii(:,:)
  integer :: np

  np = size(Lambda)
  if (np < 1) return

  allocate(Mii, source = dx/2 * this % w(1:np))
  allocate(Lii, source = 2/dx * this % L(1:np,1:np))

  call SolveGeneralizedEigenproblem(Lii, Mii, Lambda, S)

end subroutine GetEllipticEigensystem

!-------------------------------------------------------------------------------
!> Computes operators for condensed CG-SEM diffusion problem

subroutine GetEllipticSuboperators(this, dx, c, nu, Aib, Abb, Aii_inv)
  class(CG_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx           !< element length
  real(RNP), intent(in)  :: c            !< coefficient of linear term
  real(RNP), intent(in)  :: nu           !< diffusivity
  real(RNP), intent(out) :: Aib(:,:)     !< interior-boundary part, dim (po-1,2)
  real(RNP), intent(out) :: Abb(:,:)     !< boundary-boundary part, dim (2,2)
  real(RNP), intent(out) :: Aii_inv(:,:) !< Aᵢᵢ⁻¹, dimension (po-1,po-1)

  real(RNP), allocatable :: S(:,:), Lambda(:), D_inv(:)
  real(RNP) :: g0, g1
  integer   :: i, j, np

  associate(po => this%po, Ms => this%w, Ls => this%L)

    np = po - 1

    allocate(S(np,np), Lambda(np), D_inv(np))
    call this % GetEllipticEigensystem(dx, S, Lambda)

    g0 = c * dx / 2
    g1 = nu * 2 / dx

    do i = 1, np
      Aib(i,1)  =  g1 * Ls( 0,i)
      Aib(i,2)  =  g1 * Ls(po,i)
    end do

    Abb(1,1)  =  g0 * Ms( 0)  +  g1 * Ls( 0, 0)
    Abb(2,1)  =                  g1 * Ls(po, 0)
    Abb(1,2)  =                  g1 * Ls( 0,po)
    Abb(2,2)  =  g0 * Ms(po)  +  g1 * Ls(po,po)

    D_inv = 1 / (c + nu * Lambda)
    do j = 1, np
    do i = 1, np
      Aii_inv(i,j) = sum(S(i,:) * D_inv * S(j,:))
    end do
    end do

  end associate

end subroutine GetEllipticSuboperators

!===============================================================================

end module CG_Element_Operators_1D
