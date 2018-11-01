!> summary:  Element operators for continuous Galerkin-SEM
!> author:   Joerg Stiller
!> date:     2018/09/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Element operators for continuous Galerkin-SEM
!===============================================================================

module CG_Elliptic_Operator_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, ZERO
  use Execution_Control, only: Error
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Element operators for continuous Galerkin-SEM
  !>
  !> Provides the column matrix of generalized eigenvectors `S` and
  !> the diagonal matrix of eigenvalues `Λ = Lambda` such that
  !>
  !>     Sᵀ Lᵢᵢ S = Λ
  !>     Sᵀ Mᵢᵢ S = I
  !>
  !> where `Lᵢᵢ` and `Mᵢᵢ` the standard stiffness matrix and the standard
  !> diagonal mass matrix restricted to the interior points.

  type, extends(StandardOperators1D), public :: CG_EllipticOperator1D
    private
    real(RNP), allocatable :: S(:,:)    !< generalized interior eigenvectors
    real(RNP), allocatable :: Lambda(:) !< generalized interior eigenvalues
  contains
    procedure :: BuildInteriorEigensystem
    procedure :: GetInteriorEigensystem
    procedure :: GetEllipticSuboperators
    procedure :: GetStiffnessMatrix
    final     :: Delete_CG_EllipticOperator1D
  end type CG_EllipticOperator1D

contains

!-------------------------------------------------------------------------------
!> Provides the generalized eigensystem for interior stiffness and mass matrices

subroutine BuildInteriorEigensystem(this)
  class(CG_EllipticOperator1D), intent(inout) :: this

  real(RNP), allocatable :: Lii(:,:)
  integer :: np

  if (allocated(this%S)) return

  np = this%po - 1
  allocate(this % S(np,np), this % Lambda(np))

  if (np < 1) return

  allocate(Lii, source = this % L(1:np,1:np))
  associate(Mii => this % w(1:np))
    call SolveGeneralizedEigenproblem(Lii, Mii, this%Lambda, this%S)
  end associate

end subroutine BuildInteriorEigensystem

!-------------------------------------------------------------------------------
!> Returns the generalized eigenvectors and eigenvalues to the standard interior
!> stiffness and diagonal mass matrices

subroutine GetInteriorEigensystem(this, S, Lambda)
  class(CG_EllipticOperator1D), intent(inout) :: this
  real(RNP), intent(out) :: S(this%po-1,this%po-1)
  real(RNP), intent(out) :: Lambda(this%po-1)

  if (.not. allocated(this%S)) then
    call BuildInteriorEigensystem(this)
  end if

  S = this % S
  Lambda = this % Lambda

end subroutine GetInteriorEigensystem

!-------------------------------------------------------------------------------
!> Computes operators for condensed CG-SEM diffusion problem

subroutine GetEllipticSuboperators(this, dx, c, nu, Hbb, Hbi, Hii_inv)
  class(CG_EllipticOperator1D), intent(in) :: this
  real(RNP), intent(in)  :: dx                !< element length
  real(RNP), intent(in)  :: c                 !< coefficient of linear term
  real(RNP), intent(in)  :: nu                !< diffusivity
  real(RNP), intent(out) :: Hbb(2,2)          !< boundary-boundary part
  real(RNP), intent(out) :: Hbi(2,this%po-1)  !< boundary-interior part
  real(RNP), intent(out) :: Hii_inv(this%po-1,this%po-1) !< Hᵢᵢ⁻¹

  real(RNP), allocatable :: D(:)
  real(RNP) :: g0, g1
  integer   :: i, j, np

  g0 = c * dx / 2
  g1 = nu * 2 / dx

  if (.not. allocated(this%S)) then
    call Error( 'GetEllipticSuboperators'                       &
              , 'requires preceding call to BuildInteriorEigensystem' &
              , 'CG_Elliptic_Operator_1D'                       )
  end if

  associate(po => this%po, Ms => this%w, Ls => this%L, S => this%S)

    Hbb(1,1)  =  g0 * Ms( 0)  +  g1 * Ls( 0, 0)
    Hbb(2,1)  =                  g1 * Ls(po, 0)
    Hbb(1,2)  =                  g1 * Ls( 0,po)
    Hbb(2,2)  =  g0 * Ms(po)  +  g1 * Ls(po,po)

    np = po - 1

    do i = 1, np
      Hbi(1,i)  =  g1 * Ls( 0,i)
      Hbi(2,i)  =  g1 * Ls(po,i)
    end do

    allocate(D, source = 1/(g0 + g1*this%Lambda))
    do j = 1, np
    do i = 1, np
      Hii_inv(i,j) = sum(S(i,:) * S(j,:) * D)
    end do
    end do

  end associate

end subroutine GetEllipticSuboperators

!-------------------------------------------------------------------------------
!> Returns the 1D element stiffness matrix for the continuous Galerkin SEM
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:P,0:P,-1:1)`, where `P = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine GetStiffnessMatrix(this, dx, bc, Le)
  class(CG_EllipticOperator1D), intent(in) :: this
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
!> Finalization of a CG_EllipticOperator1D object

subroutine Delete_CG_EllipticOperator1D(this)
  type(CG_EllipticOperator1D), intent(inout) :: this

  if (allocated(this%S     )) deallocate(this%S     )
  if (allocated(this%Lambda)) deallocate(this%Lambda)

end subroutine Delete_CG_EllipticOperator1D

!===============================================================================

end module CG_Elliptic_Operator_1D
