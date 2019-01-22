!> summary:  Direct 1D elliptic solver for hybrid IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2019/01/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Direct 1D elliptic solver for hybrid IP/DG-SEM
!===============================================================================

module IP_Hybridized_Solver_1D
  use Kind_Parameters,  only: RNP
  use Costants,         only: ZERO
  use Linear_Equations, only: TridiagonalSolver, CyclicTridiagonalSolver
  use IP_Element_Operators_1D
  implicit none
  private

  public :: HybridEllipticSolver

contains

!-------------------------------------------------------------------------------
!> Direct elliptic solver based on hybridization

subroutine HybridEllipticSolver(eop, dx, c, nu, bc, ub, f, u, standby)
  class(IP_ElementOperators1D), intent(in) :: eop !< element operators
  real(RNP), intent(in)  :: dx       !< element width
  real(RNP), intent(in)  :: c        !< coefficient of linear term
  real(RNP), intent(in)  :: nu       !< diffusivity
  character, intent(in)  :: bc(2)    !< boundary conditions {'D','N','P'}
  real(RNP), intent(in)  :: ub(2)    !< Dirichlet boundary values
  real(RNP), intent(in)  :: f(0:,:)  !< source including Neumann BC
  real(RNP), intent(out) :: u(0:,:)  !< solution

  !> optionally keep suboperators for repeated application [F]
  logical, optional, intent(in) :: standby

  ! variables ..................................................................

  ! suboperators
  real(RNP), allocatable, save :: Aib(:,:,:)     ! interior-boundary operator Â
  real(RNP), allocatable, save :: Aii_inv(:,:,:) ! inverse interior operator Ã⁻¹

  ! flux system
  real(RNP), allocatable :: Af(:,:) ! flux system matrix
  real(RNP), allocatable :: ff(:)   ! flux RHS / solution

  real(RNP) :: dx_(-1:1), tau
  integer   :: po, ne

  ! preprocessing ............................................................

  po  = eop%po
  ne  = size(f,2)
  dx_ = dx
  tau = 2 * nu * eop%PenaltyFactor(dx)

  if (allocated(Aib)) then
    if (ubound(Aib,1) /= po) deallocate(Aib, Aii_inv)
  end if

  if (.not. allocated(Aib)) then
    allocate( Aib     (0:po, 1:2 , -1:1), source = ZERO)
    allocate( Aii_inv (0:po, 0:po, -1:1), source = ZERO)
    select case(ne)
    case(1)
      call eop % GetEllipticSuboperators( dx_, bc, c, nu,              &
                                          Aib(:,:, 0), Aii_inv(:,:, 0) )
    case default
      ! left element
      call eop % GetEllipticSuboperators( dx_, [ bc(1), ' ' ], c, nu,  &
                                          Aib(:,:,-1), Aii_inv(:,:,-1) )

      ! interior element(s)
      call eop % GetEllipticSuboperators( dx_, [ ' ', ' ' ], c, nu,    &
                                          Aib(:,:, 0), Aii_inv(:,:, 0) )
      ! right element
      call eop % GetEllipticSuboperators( dx_, [ ' ', bc(2) ], c, nu,  &
                                          Aib(:,:, 1), Aii_inv(:,:, 1) )
    end select
  end if

  ! solution ...................................................................

!#### if ne > 1:

  call BuildFluxSystem(Aib, Aii_inv, tau, bc, ub, f, Af, ff)

  ! clean-up ...................................................................

  if (present(standby)) then
    if (standby) return
  end if

  deallocate(Aib, Aii_inv)

end subroutine HybridEllipticSolver

!-------------------------------------------------------------------------------
!> Build the flux system for two or more elements

subroutine BuildFluxSystem(Aib, Aii_inv, bc, ub, f, Af, ff)
  real(RNP), intent(in) :: Aib    (0:,1:,-1:)    !< interior-boundary operators
  real(RNP), intent(in) :: Aii_inv(0:,0:,-1:)    !< inverse interior operators
  real(RNP), intent(in) :: tau                   !< penalty, τ = 2μν
  character, intent(in) :: bc(2)                 !< boundary conditions
  real(RNP), intent(in) :: ub(2)                 !< Dirichlet boundary values
  real(RNP), intent(in) :: f (0:,:)              !< source including Neumann BC
  real(RNP), allocatable, intent(out) :: Af(:,:) !< flux system matrix
  real(RNP), allocatable, intent(out) :: ff(:)   !< flux system RHS

  integer :: i, po, ne, nf

  ! intialization ..............................................................

  po = ubound(f,1)
  ne = ubound(f,2)

  if (bc(2) == 'P') then
    nf = ne
  else
    nf = ne - 1
  end if

  allocate(Af(nf,3), ff(nf))

  ! system matrix ..............................................................

  select case(nf)

  case(1) ! non-periodic, two elements

    Af(1,1) = 0
    Af(1,2) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), Aib(:,2,-1))) &
            + dot_product(Aib(:,1, 1), matmul(Aii_inv(:,:, 1), Aib(:,1, 1))) &
            - 1 / (2*tau)
    Af(1,3) = 0

  case default

    Af(1,1) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), Aib(:,1,-1)))
    Af(1,2) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), Aib(:,2,-1))) &
            + dot_product(Aib(:,1, 0), matmul(Aii_inv(:,:, 0), Aib(:,1, 0))) &
            - 1 / (2*tau)
    Af(1,3) = dot_product(Aib(:,1, 0), matmul(Aii_inv(:,:, 0), Aib(:,2, 0)))

    Af(2:nf-1,1) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,1,0)))
    Af(2:nf-1,2) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,2,0))) &
                 + dot_product(Aib(:,1,0), matmul(Aii_inv(:,:,0), Aib(:,1,0))) &
                 - 1 / (2*tau)
    Af(2:nf-1,3) = dot_product(Aib(:,1,0), matmul(Aii_inv(:,:,0), Aib(:,2,0)))

    Af(nf,1) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,1,0)))
    Af(nf,2) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,2,0))) &
             + dot_product(Aib(:,1,1), matmul(Aii_inv(:,:,1), Aib(:,1,1))) &
             - 1 / (2*tau)
    Af(nf,3) = dot_product(Aib(:,1,1), matmul(Aii_inv(:,:,1), Aib(:,2,1)))

  end select

  ! RHS ........................................................................

  select case(nf)

  case(1) ! non-periodic, two elements

    ff(1) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), f(:,1))) &
          + dot_product(Aib(:,1, 1), matmul(Aii_inv(:,:, 1), f(:,2)))

  case default

    ff(1) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), f(:,1))) &
          + dot_product(Aib(:,1, 0), matmul(Aii_inv(:,:, 0), f(:,2)))

    do i = 2, nf-1
      ff(i) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), f(:,i  ))) &
            + dot_product(Aib(:,1,0), matmul(Aii_inv(:,:,0), f(:,i+1)))
    end do

    if (nf == ne-1) then
      ! non-periodic
      ff(nf) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), f(:,ne-1))) &
             + dot_product(Aib(:,1,1), matmul(Aii_inv(:,:,1), f(:,ne  )))
    else
      ! periodic
      ff(nf) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), f(:,ne))) &
             + dot_product(Aib(:,1,1), matmul(Aii_inv(:,:,1), f(:,1 )))
    end if

end subroutine BuildFluxSystem

!===============================================================================

end module IP_Hybridized_Solver_1D
