!> summary:  Direct 1D elliptic solver for SEM with static condensation
!> author:   Joerg Stiller
!> date:     2018/09/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Direct 1D elliptic solver for SEM with static condensation
!===============================================================================

module CG_Condensed_Solver_1D
  use Kind_Parameters,  only: RNP
  use Linear_Equations, only: TridiagonalSolver, CyclicTridiagonalSolver
  use CG_Elliptic_Operator_1D
  implicit none
  private

  public :: CondensedEllipticSolver

contains

!-------------------------------------------------------------------------------
!> Direct elliptic solver based on static condensation
!>
!> Element suboperators provided with `eop` must be initialized by call to
!> prior call to `eop % BuildInteriorEigensystem`.

subroutine CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
  class(CG_EllipticOperator1D), intent(in) :: eop !< element operators
  real(RNP), intent(in)    :: dx       !< element width
  real(RNP), intent(in)    :: c        !< coefficient of linear term
  real(RNP), intent(in)    :: nu       !< diffusivity
  character, intent(in)    :: bc(:)    !< boundary conditions {'D','N','P'}
  real(RNP), intent(in)    :: f(0:,:)  !< source including Neumann BC
  real(RNP), intent(inout) :: u(0:,:)  !< solution with Dirichlet BC set

  ! variables ..................................................................

  ! suboperators
  real(RNP), allocatable :: Hbb(:,:)     ! boundary-boundary element operator
  real(RNP), allocatable :: Hbi(:,:)     ! boundary-interior element operator
  real(RNP), allocatable :: Hii_inv(:,:) ! inverse of interior-interior part

  ! condensed system
  real(RNP), allocatable :: Ac(:,:)      ! condensed system matrix
  real(RNP), allocatable :: fc(:)        ! condensed RHS / solution

  ! preprocessing ..............................................................

  ! suboperators required for condensation
  allocate(Hbb(2,2), Hbi(2,eop%po-1), Hii_inv(eop%po-1,eop%po-1))
  call eop % GetEllipticSuboperators(dx, c, nu, Hbb, Hbi, Hii_inv)

  ! static condensation
  call BuildCondensedSystem(Hbb, Hbi, Hii_inv, bc, f, u, Ac, fc)

  ! solution ...................................................................

  call SolveCondensedSystem(Ac, fc)
  call SolveElementSystems(Hbi, Hii_inv, bc, fc, f, u)

end subroutine CondensedEllipticSolver

!-------------------------------------------------------------------------------
!> Build the condensed system

subroutine BuildCondensedSystem(Hbb, Hbi, Hii_inv, bc, f, u, Ac, fc)
  real(RNP), intent(in) :: Hbb(:,:)     !< boundary-boundary op.
  real(RNP), intent(in) :: Hbi(:,:)     !< boundary-interior op.
  real(RNP), intent(in) :: Hii_inv(:,:) !< inverse interior op.
  character, intent(in) :: bc(2)        !< boundary conditions
  real(RNP), intent(in) :: f(0:,:)      !< source including Neumann BC
  real(RNP), intent(in) :: u(0:,:)      !< initial values including Dirichlet BC
  real(RNP), allocatable, intent(out) :: Ac(:,:) !< condensed system matrix
  real(RNP), allocatable, intent(out) :: fc(:)   !< condensed RHS

  ! variables ..................................................................

  real(RNP) :: a(2,2), b(2,size(Hbi,2))
  real(RNP), allocatable :: fci(:,:)
  integer   :: i, i1, i2, po, ne, np

  ! intialization ..............................................................

  po = ubound(f,1)
  ne = ubound(f,2)
  np = po - 1

  ! left boundary
  select case(bc(1))
  case('D','P') ! exclude boundary in case of Dirichlet or periodic BC
    i1 = 1
  case default  ! include boundary in case of Neumann BC
    i1 = 0
  end select

  ! right boundary
  select case(bc(2))
  case('D')     ! exclude boundary in case of Dirichlet BC
    i2 = ne-1
  case default  ! include boundary in case of Neumann or periodic BC
    i2 = ne
  end select

  ! condensed element operators ................................................

  b = -matmul(Hbi, Hii_inv)
  a =  Hbb + matmul(b, transpose(Hbi))

  ! condensed system matrix ....................................................

  allocate(Ac(i1:i2,3))

  Ac(:,1) = a(2,1)
  Ac(:,2) = a(1,1) + a(2,2)
  Ac(:,3) = a(1,2)

  ! left boundary
  select case(bc(1))
  case('D')
    Ac(i1,1) = 0
  case('N')
    Ac(i1,1) = 0
    Ac(i1,2) = Ac(i1,2) - a(2,2)
  end select

  ! right boundary
  select case(bc(2))
  case('D')
    Ac(i2,3) = 0
  case('N')
    Ac(i2,2) = Ac(i2,2) - a(1,1)
    Ac(i2,3) = 0
  end select

  ! condensed RHS ..............................................................

  allocate(fc(i1:i2))

  ! extract element boundary values
  if (bc(1) == 'N') then
    fc( 0) = f( 0, 1)
  end if
  do i = 1, i2
    fc(i) = f(po,i)
  end do

  ! compute interior contributions
  fci = matmul(b, f(1:po-1,:))

  ! add interior contributions at left boundary
  select case(bc(1))
  case('N')
    fc(0) = fc(0) + fci(1,1)
  end select

  ! add interior contributions at inner points
  do i = 1, ne-1
    fc(i) = fc(i) + fci(2,i) + fci(1,i+1)
  end do

  ! add interior contributions at right boundary
  select case(bc(2))
  case('N')
    fc(ne) = fc(ne) + fci(2,ne)
  case('P')
    fc(ne) = fc(ne) + fci(2,ne) + fci(1,1)
  end select

  ! Dirichlet contributions
  if (bc(1) == 'D') fc(i1) = fc(i1) - a(2,1) * u( 0, 1)
  if (bc(2) == 'D') fc(i2) = fc(i2) - a(1,2) * u(po,ne)

end subroutine BuildCondensedSystem

!------------------------------------------------------------------------------
!> Solution of the condensed system using tridiagonal Gauss elimination

subroutine SolveCondensedSystem(Ac, fc)
  real(RNP), intent(inout) :: Ac(:,:) !< system matrix, destroyed on output
  real(RNP), intent(inout) :: fc(:)   !< RHS (in) / solution (out)

  logical :: regular
  integer :: n

  n = size(fc)

  if (n == 0) then
    return

  else if (n == 1) then
    fc(1) = fc(1) / Ac(1,2)

  else
    associate(a => Ac(:,1), b => Ac(:,2), c => Ac(:,3))

      regular = abs(b(1)) > abs(c(1)) .or. abs(b(n)) > abs(a(n))
      if (a(1) == 0 .and. c(n) == 0) then
        call TridiagonalSolver(fc, a, b, c, regular)
      else
        call CyclicTridiagonalSolver(fc, a, b, c, regular)
      end if

    end associate

  end if

end subroutine SolveCondensedSystem

!------------------------------------------------------------------------------
!> Expand condensed solution by solution of element systems

subroutine SolveElementSystems(Hbi, Hii_inv, bc, uc, f, u)
  real(RNP), intent(in)    :: Hbi(:,:)     !< boundary-interior element op.
  real(RNP), intent(in)    :: Hii_inv(:,:) !< inverse interior element op.
  character, intent(in)    :: bc(2)        !< boundary conditions
  real(RNP), intent(in)    :: uc(:)        !< condensed solution
  real(RNP), intent(in)    :: f(0:,:)      !< full RHS
  real(RNP), intent(inout) :: u(0:,:)      !< full solution

  integer :: i, k, po, ne

  ! intialization ..............................................................

  po = ubound(u,1)
  ne = ubound(u,2)

  ! inject condensed solution ..................................................

  ! left boundary
  select case(bc(1))
  case('N')
    u(0,1) = uc(1)
  case('P')
    u(0,1) = uc(size(uc))
  end select

  ! start index
  select case(bc(1))
  case('N') ! Neumann BC: skip uc(k=1) since already injected to left boundary
    k = 2
  case default
    k = 1
  end select

  ! interior
  do i = 1, ne-1
    u(po, i  ) = uc(k)
    u( 0, i+1) = uc(k)
    k = k + 1
  end do

  ! right boundary
  select case(bc(2))
  case('N','P')
    u(po,ne) = uc(size(uc))
  end select

  ! solve interior subsystems ..................................................

  k = po-1
  do i = 1, ne
    u(1:k,i) = matmul(Hii_inv, f(1:k,i) - (Hbi(1,:)*u(0,i) + Hbi(2,:)*u(po,i)))
  end do

end subroutine SolveElementSystems

!===============================================================================

end module CG_Condensed_Solver_1D
