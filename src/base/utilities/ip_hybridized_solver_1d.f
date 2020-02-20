!> summary:  Direct 1D elliptic solver for hybrid IP/DG-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2019/01/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Direct 1D elliptic solver for hybrid IP/DG-SEM
!===============================================================================

module IP_Hybridized_Solver_1D
  use Kind_Parameters,  only: RNP
  use Constants,        only: ZERO
  use Linear_Equations, only: TridiagonalSolver, CyclicTridiagonalSolver
  use IP_Element_Operators_1D
  implicit none
  private

  public :: HybridEllipticSolver

  interface HybridEllipticSolver
    procedure :: HybridEllipticSolver__w_svv
    procedure :: HybridEllipticSolver__n_svv
  end interface

contains

!-------------------------------------------------------------------------------
!> Direct elliptic solver based on hybridization including SVV

subroutine HybridEllipticSolver__w_svv( eop, dx, c, nu, nu_svv &
                                      , bc, f, u, standby      )
  class(IP_ElementOperators1D), intent(in) :: eop !< element operators
  real(RNP), intent(in)  :: dx       !< element width
  real(RNP), intent(in)  :: c        !< coefficient of linear term
  real(RNP), intent(in)  :: nu       !< diffusivity
  real(RNP), intent(in)  :: nu_svv   !< SVV diffusivity [0]
  character, intent(in)  :: bc(2)    !< boundary conditions {'D','N','P'}
  real(RNP), intent(in)  :: f(0:,:)  !< source including Neumann BC
  real(RNP), intent(out) :: u(0:,:)  !< solution

  !> optionally keep suboperators for repeated application [F]
  logical, optional, intent(in) :: standby

  ! variables ..................................................................

  ! suboperators
  real(RNP), allocatable, save :: Aib(:,:,:)     ! interior-boundary operator Â
  real(RNP), allocatable, save :: Aii_inv(:,:,:) ! inverse interior operator Ã⁺

  ! flux system
  real(RNP), allocatable :: Af(:,:) ! flux system matrix
  real(RNP), allocatable :: ff(:)   ! flux RHS / solution

  real(RNP) :: dx_(-1:1), tau
  integer   :: po, ne

  ! preprocessing ............................................................

  po  = eop%po
  ne  = size(f,2)
  dx_ = dx

  tau = 2 * (nu + nu_svv) * eop%PenaltyFactor(dx)

  if (allocated(Aib)) then
    if (ubound(Aib,1) /= po) deallocate(Aib, Aii_inv)
  end if

  if (.not. allocated(Aib)) then
    allocate( Aib     (0:po, 1:2 , -1:1), source = ZERO)
    allocate( Aii_inv (0:po, 0:po, -1:1), source = ZERO)
    select case(ne)
    case(1)
      call eop % Get_EllipticSuboperators( dx_, bc,             c, nu, nu_svv, &
                                           Aib(:,:, 0), Aii_inv(:,:, 0)        )
    case default
      ! left element
      call eop % Get_EllipticSuboperators( dx_, [ bc(1), ' ' ], c, nu, nu_svv, &
                                           Aib(:,:,-1), Aii_inv(:,:,-1)        )

      ! interior element(s)
      call eop % Get_EllipticSuboperators( dx_, [ ' ', ' ' ],   c, nu, nu_svv, &
                                           Aib(:,:, 0), Aii_inv(:,:, 0)        )
      ! right element
      call eop % Get_EllipticSuboperators( dx_, [ ' ', bc(2) ], c, nu, nu_svv, &
                                           Aib(:,:, 1), Aii_inv(:,:, 1)        )
    end select
  end if

  ! solution ...................................................................

  if (ne > 1) then
    call BuildFluxSystem(Aib, Aii_inv, tau, bc, f, Af, ff)
    call SolveFluxSystem(Af, ff)
    call SolveElementSystems(Aib, Aii_inv, bc, ff, f, u)
  else
    u = matmul(Aii_inv(:,:,0), f)
  end if

  ! clean-up ...................................................................

  if (present(standby)) then
    if (standby) return
  end if

  deallocate(Aib, Aii_inv)

end subroutine HybridEllipticSolver__w_svv

!-------------------------------------------------------------------------------
!> Direct elliptic solver based on hybridization without SVV

subroutine HybridEllipticSolver__n_svv(eop, dx, c, nu, bc, f, u, standby)
  class(IP_ElementOperators1D), intent(in) :: eop !< element operators
  real(RNP), intent(in)  :: dx       !< element width
  real(RNP), intent(in)  :: c        !< coefficient of linear term
  real(RNP), intent(in)  :: nu       !< diffusivity
  character, intent(in)  :: bc(2)    !< boundary conditions {'D','N','P'}
  real(RNP), intent(in)  :: f(0:,:)  !< source including Neumann BC
  real(RNP), intent(out) :: u(0:,:)  !< solution

  !> optionally keep suboperators for repeated application [F]
  logical, optional, intent(in) :: standby

  call HybridEllipticSolver__w_svv(eop, dx, c, nu, ZERO, bc, f, u, standby)

end subroutine HybridEllipticSolver__n_svv

!-------------------------------------------------------------------------------
!> Build the flux system for two or more elements

subroutine BuildFluxSystem(Aib, Aii_inv, tau, bc, f, Af, ff)
  real(RNP), intent(in) :: Aib    (0:,1:,-1:)    !< interior-boundary operators
  real(RNP), intent(in) :: Aii_inv(0:,0:,-1:)    !< inverse interior operators
  real(RNP), intent(in) :: tau                   !< penalty, τ = 2μν
  character, intent(in) :: bc(2)                 !< boundary conditions
  real(RNP), intent(in) :: f (0:,:)              !< source including Neumann BC
  real(RNP), allocatable, intent(out) :: Af(:,:) !< flux system matrix
  real(RNP), allocatable, intent(out) :: ff(:)   !< flux system RHS

  integer :: i, ne, nf

  ! intialization ..............................................................

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
            - 2*tau
    Af(1,3) = 0

  case default

    Af(1,1) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), Aib(:,1,-1)))
    Af(1,2) = dot_product(Aib(:,2,-1), matmul(Aii_inv(:,:,-1), Aib(:,2,-1))) &
            + dot_product(Aib(:,1, 0), matmul(Aii_inv(:,:, 0), Aib(:,1, 0))) &
            - 2*tau
    Af(1,3) = dot_product(Aib(:,1, 0), matmul(Aii_inv(:,:, 0), Aib(:,2, 0)))

    Af(2:nf-1,1) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,1,0)))
    Af(2:nf-1,2) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,2,0))) &
                 + dot_product(Aib(:,1,0), matmul(Aii_inv(:,:,0), Aib(:,1,0))) &
                 - 2*tau
    Af(2:nf-1,3) = dot_product(Aib(:,1,0), matmul(Aii_inv(:,:,0), Aib(:,2,0)))

    Af(nf,1) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,1,0)))
    Af(nf,2) = dot_product(Aib(:,2,0), matmul(Aii_inv(:,:,0), Aib(:,2,0))) &
             + dot_product(Aib(:,1,1), matmul(Aii_inv(:,:,1), Aib(:,1,1))) &
             - 2*tau
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

  end select

end subroutine BuildFluxSystem

!------------------------------------------------------------------------------
!> Solution of the condensed system using tridiagonal Gauss elimination

subroutine SolveFluxSystem(Af, ff)
  real(RNP), intent(inout) :: Af(:,:) !< system matrix, destroyed on output
  real(RNP), intent(inout) :: ff(:)   !< RHS (in) / solution (out)

  logical :: regular
  integer :: n

  n = size(ff)

  if (n == 0) then
    return

  else if (n == 1) then
    if (abs(Af(1,2)) > epsilon(Af)) then
      ff(1) = ff(1) / Af(1,2)
    else
      ff(1) = 0
    end if

  else
    associate(a => Af(:,1), b => Af(:,2), c => Af(:,3))

      regular = abs(b(1)) > abs(c(1)) .or. abs(b(n)) > abs(a(n))
      if (a(1) == 0 .and. c(n) == 0) then
        call TridiagonalSolver(ff, a, b, c, regular)
      else
        call CyclicTridiagonalSolver(ff, a, b, c, regular)
      end if

    end associate

  end if

end subroutine SolveFluxSystem

!------------------------------------------------------------------------------
!> Solution of the element systems

subroutine SolveElementSystems(Aib, Aii_inv, bc, uf, f, u)
  real(RNP), intent(in)  :: Aib    (0:,1:,-1:) !< interior-boundary operators
  real(RNP), intent(in)  :: Aii_inv(0:,0:,-1:) !< inverse interior operators
  character, intent(in)  :: bc(2)              !< boundary conditions
  real(RNP), intent(in)  :: uf(:)              !< fluxes
  real(RNP), intent(in)  :: f(0:,:)            !< full RHS
  real(RNP), intent(out) :: u(0:,:)            !< full solution

  integer :: i, ne

  ne = ubound(u, 2)

  ! left
  if (bc(1) == 'P') then
    u(:,1) = matmul( Aii_inv(:,:, 0), f(:,1) - Aib(:,1, 0) * uf(ne) &
                                             - Aib(:,2, 0) * uf(1)  )
  else
    u(:,1) = matmul( Aii_inv(:,:,-1), f(:,1) - Aib(:,2,-1) * uf(1) )
  end if

  ! interior
  do i = 2, ne-1
    u(:,i) = matmul( Aii_inv(:,:, 0), f(:,i) - Aib(:,1, 0) * uf(i-1) &
                                             - Aib(:,2, 0) * uf(i)   )
  end do

  ! right
  if (bc(2) == 'P') then
    u(:,ne) = matmul( Aii_inv(:,:, 0), f(:,ne) - Aib(:,1, 0) * uf(ne-1) &
                                               - Aib(:,2, 0) * uf(ne)   )
  else
    u(:,ne) = matmul( Aii_inv(:,:, 1), f(:,ne) - Aib(:,1, 1) * uf(ne-1) )
  end if

end subroutine SolveElementSystems

!===============================================================================

end module IP_Hybridized_Solver_1D
