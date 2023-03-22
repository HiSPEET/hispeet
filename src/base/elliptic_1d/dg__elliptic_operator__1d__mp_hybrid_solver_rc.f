!> summary:  Direct 1D elliptic solver for hybrid IP/DG-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2019/01/20, revised 2022/03/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__1D) MP_Hybrid_Solver_RC
  use Linear_Equations, only: TridiagonalSolver, CyclicTridiagonalSolver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Direct elliptic solver based on hybridization including SVV

  module subroutine HybridSolver_RC(this, dx, lambda, nu, f, bv, u, standby)

    class(DG_EllipticOperator_1D), intent(in)  :: this
    real(RNP),             intent(in)    :: dx       !< element width
    real(RNP),             intent(in)    :: lambda   !< λ
    real(RNP),             intent(in)    :: nu       !< diffusivity, ν = νᵖ+νˢ
    real(RNP),             intent(in)    :: bv(2)    !< boundary values, u or u'
    real(RNP), contiguous, intent(inout) :: f(0:,:)  !< source, being destroyed
    real(RNP), contiguous, intent(out)   :: u(0:,:)  !< solution

    !> optionally keep suboperators for repeated application [F]
    logical, optional, intent(in) :: standby

    ! variables ................................................................

    ! suboperators
    real(RNP), allocatable, save :: Aib(:,:,:)     ! interior-boundary operator Â
    real(RNP), allocatable, save :: Aii_inv(:,:,:) ! inverse interior operator Ã⁺

    ! flux system
    real(RNP), allocatable :: Af(:,:) ! flux system matrix
    real(RNP), allocatable :: ff(:)   ! flux RHS / solution

    real(RNP) :: nu_p, nu_s
    real(RNP) :: dx_(-1:1), tau
    character :: bc_sub(2)
    integer   :: po, ne

    associate(eop => this % eop, bc => this % bc)

      ! preprocessing ..........................................................

      po  = eop%po
      ne  = size(f,2)
      dx_ = dx

      nu_p = this % PhysicalDiffusivity(nu)
      nu_s = this % SpectralDiffusivity(nu)

      tau = 2 * (nu_p + nu_s) * eop%PenaltyFactor(dx)

      !$omp master

      if (allocated(Aib)) then
        if (ubound(Aib,1) /= po) deallocate(Aib, Aii_inv)
      end if

      if (.not. allocated(Aib)) then
        allocate( Aib     (0:po, 1:2 , -1:1), source = ZERO)
        allocate( Aii_inv (0:po, 0:po, -1:1), source = ZERO)
        select case(ne)
        case(1)
          call eop % Get_EllipticSuboperators( dx_, bc, lambda, nu_p, nu_s  &
                                             , Aib(:,:, 0), Aii_inv(:,:, 0) )
        case default
          ! left element
          bc_sub = [ bc(1), ' ' ]
          call eop % Get_EllipticSuboperators( dx_, bc_sub, lambda, nu_p, nu_s &
                                             , Aib(:,:,-1), Aii_inv(:,:,-1)    )

          ! interior element(s)
          bc_sub = [ ' ', ' ' ]
          call eop % Get_EllipticSuboperators( dx_, bc_sub, lambda, nu_p, nu_s &
                                             , Aib(:,:, 0), Aii_inv(:,:, 0)    )
          ! right element
          bc_sub = [ ' ', bc(2) ]
          call eop % Get_EllipticSuboperators( dx_, bc_sub, lambda, nu_p, nu_s &
                                             , Aib(:,:, 1), Aii_inv(:,:, 1)    )
        end select
      end if

      call ApplyBoundaryConditions(eop, dx, bc, bv, f)
      if (lambda == 0 .and. (all(bc == 'N') .or. all(bc == 'P'))) then
        f = f - sum(f) / size(f)
      end if

      !$omp end master
      !$omp barrier

      ! solution ...............................................................

      if (ne > 1) then
        call BuildFluxSystem(Aib, Aii_inv, tau, bc, f, Af, ff)
        call SolveFluxSystem(Af, ff)
        call SolveElementSystems(Aib, Aii_inv, bc, ff, f, u)
      else
        u = matmul(Aii_inv(:,:,0), f)
      end if

      ! clean-up ...............................................................

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(Aib, Aii_inv)
      !$omp end master

    end associate

  end subroutine HybridSolver_RC

  !-----------------------------------------------------------------------------
  !> Apply boundary conditions to RHS

  subroutine ApplyBoundaryConditions(eop, dx, bc, bv, f)
    class(DG_ElementOperators_1D), intent(in) :: eop !< IP-H element operators
    real(RNP), intent(in)    :: dx      !< element extension
    character, intent(in)    :: bc(2)   !< boundary conditions
    real(RNP), intent(in)    :: bv(2)   !< boundary values
    real(RNP), intent(inout) :: f(0:,:) !< RHS

    real(RNP), allocatable :: delta_0(:), delta_P(:)
    real(RNP) :: tau
    integer   :: ne

    ne = ubound(f,2)

    associate(po => eop%po, Ms => eop%w, Ds => eop%D)

      tau = 2 * eop % PenaltyFactor(dx)

      allocate(delta_0(0:po), source = ZERO)
      delta_0(0) = ONE

      allocate(delta_P(0:po), source = ZERO)
      delta_P(po) = ONE

      ! left boundary
      if (bc(1) == 'D') then
        f(:,1) = f(:,1) + (2/dx * Ds(0,:) + tau * delta_0) * bv(1)
      else if (bc(1) == 'N') then
        f(0,1) = f(0,1) - bv(1)
      end if

      ! right boundary
      if (bc(2) == 'D') then
        f(:,ne) = f(:,ne) + (-2/dx * Ds(po,:) + tau * delta_P) * bv(2)
      else if (bc(2) == 'N') then
        f(po,ne) = f(po,ne) + bv(2)
      end if

    end associate

  end subroutine ApplyBoundaryConditions

  !---------------------------------------------------------------------------
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

    ! intialization ............................................................

    ne = ubound(f,2)

    if (bc(2) == 'P') then
      nf = ne
    else
      nf = ne - 1
    end if

    allocate(Af(nf,3), ff(nf))

    ! system matrix ............................................................

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

    ! RHS ......................................................................

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

  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
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

  !=============================================================================

end submodule MP_Hybrid_Solver_RC
