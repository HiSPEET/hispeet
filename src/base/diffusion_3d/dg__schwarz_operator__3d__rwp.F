  public :: DG_SchwarzOperator_3D_RWP

  !-----------------------------------------------------------------------------
  !> Schwarz operator

  type, extends(DG_SchwarzOperator_3D) :: DG_SchwarzOperator_3D_RWP

    real(RWP), allocatable :: S(:,:,:)       !< eigenvectors per config
    real(RWP), allocatable :: V(:,:)         !< eigenvalues  per config
    real(RWP), allocatable :: W(:,:)         !< weights      per config
    real(RWP), allocatable :: D_inv(:,:,:,:) !< inverse 3D eigenvalues

  contains

    procedure :: InitOperators
    procedure :: InitSubdomains_CI
    procedure :: InitSubdomains_CS
    procedure :: InitSubdomains_VI

  end type DG_SchwarzOperator_3D_RWP

  ! constructor interface
  interface DG_SchwarzOperator_3D_RWP
    module procedure New_CI
    module procedure New_CS
    module procedure New_VI
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor with constant physical diffusivity

  function New_CI(eop, no, weighting, mesh, lambda, nu, bc) result(this)

    class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
    integer,        intent(in) :: no        !< number of overlapped points
    integer,        intent(in) :: weighting !< weighting method {0,1,3,5,7,9}
    class(Mesh_3D), intent(in) :: mesh      !< mesh partition
    real(RNP),      intent(in) :: lambda    !< Helmholtz parameter
    real(RNP),      intent(in) :: nu        !< physical diffusivity
    character,      intent(in) :: bc(:)     !< BC {'D','N','P'}

    type(DG_SchwarzOperator_3D_RWP) :: this

    call this % InitOperators(eop, no, weighting)
    call this % InitSubdomains(mesh, lambda, nu, bc)

  end function New_CI

  !-----------------------------------------------------------------------------
  !> Constructor with constant physical and spectral diffusivities

  function New_CS(eop, no, weighting, mesh, lambda, nu, nu_svv, bc) result(this)

    class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
    integer,        intent(in) :: no        !< number of overlapped points
    integer,        intent(in) :: weighting !< weighting method {0,1,3,5,7,9}
    class(Mesh_3D), intent(in) :: mesh      !< mesh partition
    real(RNP),      intent(in) :: lambda    !< Helmholtz parameter
    real(RNP),      intent(in) :: nu        !< physical diffusivity
    real(RNP),      intent(in) :: nu_svv    !< spectral diffusivity
    character,      intent(in) :: bc(:)     !< BC {'D','N','P'}

    type(DG_SchwarzOperator_3D_RWP) :: this

    real(RNP) :: svv

    svv = nu_svv / (nu + nu_svv)

    call this % InitOperators(eop, no, weighting, svv)
    call this % InitSubdomains(mesh, lambda, nu, nu_svv, bc)

  end function New_CS

  !-----------------------------------------------------------------------------
  !> Constructor with variable diffusivity

  function New_VI(eop, no, weighting, mesh, lambda, nu, bc) result(this)

    class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
    integer,        intent(in) :: no          !< number of overlapped points
    integer,        intent(in) :: weighting   !< weighting method {0,1,3,5,7,9}
    class(Mesh_3D), intent(in) :: mesh        !< mesh partition
    real(RNP),      intent(in) :: lambda      !< Helmholtz parameter
    real(RNP),      intent(in) :: nu(:,:,:,:) !< diffusivity
    character,      intent(in) :: bc(:)       !< BC {'D','N','P'}

    type(DG_SchwarzOperator_3D_RWP) :: this

    call this % InitOperators(eop, no, weighting)
    call this % InitSubdomains(eop, mesh, lambda, nu, bc)

  end function New_VI

  !-----------------------------------------------------------------------------
  !> Build the 1D eigenvalues, eigenvectors and weights

  subroutine InitOperators(this, eop, no, weighting, svv)
    class(DG_SchwarzOperator_3D_RWP), intent(inout) :: this
    class(DG_ElementOperators_1D), intent(in) :: eop !< DG-SE operators
    integer,   intent(in) :: no        !< number of overlapped points
    integer,   intent(in) :: weighting !< weighting method {0,1,3,5,7,9}
    real(RNP), intent(in), optional :: svv !< ratio νˢ/(ν + νˢ)

    ! local variables ..........................................................

    real(RNP) :: Ws(eop%po + 1 + max(2*no,0))
    character :: bc(2)

    integer :: nb, nc, ns
    integer :: i, j, k

    ! initialization ...........................................................

    nb = size(DG_SCHWARZ_BC_3D) ! number of supported boundary conditions
    nc = nb ** 2                ! number of 1D boundary configurations
    ns = size(Ws)               ! number of subdomain points per direction

    !$omp single

    this % po = eop % po
    this % no = max(0,no)
    this % nc = nc
    this % restrictive = weighting == 9 ! top-hat weighting !

    if (allocated(this % S)) deallocate(this % S)
    if (allocated(this % V)) deallocate(this % V)
    if (allocated(this % W)) deallocate(this % W)

    allocate(this % S(ns,ns,nc), source = ZERO)
    allocate(this % V(ns,nc)   , source = ONE )
    allocate(this % W(ns,nc)   , source = ZERO)

    !$omp end single

    ! eigensystems and weights .................................................

    call WeightDistribution(eop%x, this%no, weighting, Ws)

    !$omp do collapse(2)
    do j = 1, nb
    do i = 1, nb

      k = i + nb * (j - 1)

      bc(1) = DG_SCHWARZ_BC_3D(i)
      bc(2) = DG_SCHWARZ_BC_3D(j)

      call ComputeGenericOperators( eop, svv, this%no, bc, Ws &
                                  , S  = this % S(:,:,k)      &
                                  , V  = this % V(:,k)        &
                                  , W  = this % W(:,k)        )
    end do
    end do

  end subroutine InitOperators

  !-----------------------------------------------------------------------------
  !> Computes the eigenvectors, eigenvalues and weights for a 1D subdomain
  !>
  !> This procedure provides the eigenvectors and eigenvalues for 1D subdomains
  !> assuming a constant element width of dx=1. The actual element width is
  !> considered be adjusting the coefficients c0, c1, c2, and c3, as detailed in
  !> the description of `DG_SchwarzOperator_3D`.
  !> The argument `bc` specifies the left and right boundary conditions.
  !> `D`indicates a Dirichlet boundary and `N` a Neumann boundary.
  !> Otherwise, the existence of a neighbor element is assumed.
  !>
  !> Although no exterior node layers exist at boundaries, the corresponding
  !> entries are retained for regularity and set to `0`.

  subroutine ComputeGenericOperators(eop, svv, no, bc, Ws, S, V, W)
    class(DG_ElementOperators_1D), intent(in) :: eop !< IP-DG element operators
    real(RNP),  intent(in), optional :: svv !< ratio νˢ/(ν + νˢ) [0]
    integer,    intent(in)  :: no           !< overlap
    character,  intent(in)  :: bc(2)        !< left/right boundary conditions
    real(RNP),  intent(in)  :: Ws(-no:)     !< standard weights
    real(RNP),  intent(out) :: S(-no:,-no:) !< subdomain eigenvectors
    real(RNP),  intent(out) :: V(-no:)      !< subdomain eigenvalues
    real(RNP),  intent(out) :: W(-no:)      !< subdomain weights

    character, parameter   :: ii(2) = [ ' ', ' ' ]
    real(RNP), parameter   :: dx(-1:1) = 1
    real(RNP), allocatable :: Le_ii(:,:,:), Le_bc(:,:,:)
    real(RDP), allocatable :: Ld(:,:), Md(:), Sd(:,:), vd(:)

    real(RNP) :: nu, nu_svv
    integer   :: po, np
    integer   :: k

    po = eop%po
    np = po + 1

    allocate(Le_ii(0:po,0:po,-1:1))
    allocate(Le_bc(0:po,0:po,-1:1))

    if (present(svv)) then
      nu_svv = svv
    else
      nu_svv = ZERO
    end if
    nu = max(ONE - nu_svv, ZERO)

    ! stiffness matrix for interior element with dx=1, bc = ii
    call eop % Get_DiffusionMatrix(dx, ii, nu, nu_svv, Ae = Le_ii)

    ! stiffness matrix for given boundary conditions
    call eop % Get_DiffusionMatrix(dx, bc, nu, nu_svv, Ae = Le_bc)

    ! initialization of eigenvectors S, eigenvalues V and weights W
    S = 0
    V = 0
    W = 0

    associate(Ms => eop % W)

      if (all(bc == ' ')) then

        !-----------------------------------------------------------------------
        ! interior-interior configuration

        ! mass matrix ..........................................................

        allocate(Md(-no:po+no))

        Md(-no:-1   ) = real( HALF * Ms( np-no:po   ), kind = RDP )
        Md(  0:po   ) = real( HALF * Ms(     0:po   ), kind = RDP )
        Md( np:po+no) = real( HALF * Ms(     0:no-1 ), kind = RDP )

        ! stiffness matrix .....................................................

        allocate(Ld(-no:po+no,-no:po+no), source = ZERO)

        ! lines from preceding element
        Ld(-no:-1, -no:-1) = real( Le_ii(np-no:po, np-no:po,  0), RDP )
        Ld(-no:-1,   0:po) = real( Le_ii(np-no:po,     0:po,  1), RDP )

        ! lines from present element
        Ld(0:po, -no:-1   ) = real(Le_ii(0:po, np-no:po  , -1), RDP )
        Ld(0:po,   0:po   ) = real(Le_ii(0:po,     0:po  ,  0), RDP )
        Ld(0:po,  np:po+no) = real(Le_ii(0:po,     0:no-1,  1), RDP )

        ! lines from following element
        Ld(np:po+no,  0:po   ) = real(Le_ii(0:no-1, 0:po  , -1), RDP )
        Ld(np:po+no, np:po+no) = real(Le_ii(0:no-1, 0:no-1,  0), RDP )

        ! eigenvectors and eigenvalues .........................................

        allocate(Sd, mold = Ld)
        allocate(Vd, mold = Md)

        call SolveGeneralizedEigenproblem(Ld, Md, Vd, Sd)

        ! inject eigenvectors and eigenvalues
        S = real( Sd, kind = RWP )
        V = real( Vd, kind = RWP )

        ! weights ..............................................................

        W = real( Ws, kind = RWP )

      else if (all(bc /= ' ')) then

        !-----------------------------------------------------------------------
        ! boundary-boundary configuration

        ! mass matrix ..........................................................

        allocate(Md(0:po))

        Md(0:po) = real( HALF * Ms, kind = RDP )

        ! stiffness matrix .....................................................

        allocate(Ld(0:po,0:po))

        Ld(0:po,0:po) = real( Le_bc(:,:,0), kind = RDP )

        ! eigenvectors and eigenvalues .........................................

        allocate(Sd(0:po,0:po), Vd(0:po))

        call SolveGeneralizedEigenproblem(Ld, Md, Vd, Sd)

        ! inject eigenvectors and eigenvalues
        S(0:po,0:po) = real( Sd, kind = RWP )
        V(0:po)      = real( Vd, kind = RWP )

        ! weights ..............................................................

        W(0:po) = real( Ws(0:po), kind = RWP )

        ! left boundary zone
        do k = 0, no-1
          W(k) = W(k) + real(Ws(-k-1), RWP)
        end do

        ! right boundary zone
        do k = po-no+1, po
          W(k) = W(k) + real(Ws(2*po + 1 - k), RWP)
        end do

      else if (bc(1) /= ' ') then

        !-----------------------------------------------------------------------
        ! boundary-interior configuration

        ! mass matrix ..........................................................

        allocate(Md(0:po+no))

        Md(  0:po   ) = real( HALF * Ms(     0:po   ), kind = RDP )
        Md( np:po+no) = real( HALF * Ms(     0:no-1 ), kind = RDP )

        ! stiffness matrix .....................................................

        allocate(Ld(0:po+no,0:po+no), source = 0D0)

        ! lines from present element
        Ld( 0:po,  0:po   ) = real( Le_bc(0:po, 0:po  ,  0), RDP )
        Ld( 0:po, np:po+no) = real( Le_bc(0:po, 0:no-1,  1), RDP )

        ! lines from following element
        Ld(np:po+no,  0:po   ) = real(Le_ii(0:no-1, 0:po  , -1), RDP )
        Ld(np:po+no, np:po+no) = real(Le_ii(0:no-1, 0:no-1,  0), RDP )

        ! eigenvectors and eigenvalues .........................................

        allocate(Sd(0:po+no,0:po+no), Vd(0:po+no))

        call SolveGeneralizedEigenproblem(Ld, Md, Vd, Sd)

        ! inject eigenvectors and eigenvalues

        S(0:po+no,0:po+no) = real( Sd, kind = RWP )
        V(0:po+no)         = real( Vd, kind = RWP )

        ! weights ..............................................................

        W(0:) = real( Ws(0:), kind = RWP )

        ! boundary zone
        do k = 0, no-1
          W(k) = W(k) + real(Ws(-k-1), RWP)
        end do

      else

        !-----------------------------------------------------------------------
        ! interior-boundary configuration

        ! mass matrix ..........................................................

        allocate(Md(-no:po))

        Md(-no:-1   ) = real( HALF * Ms( np-no:po   ), kind = RDP )
        Md(  0:po   ) = real( HALF * Ms(     0:po   ), kind = RDP )

        ! stiffness matrix .....................................................

        allocate(Ld(-no:po,-no:po), source = 0D0)

        ! lines from preceding element
        Ld(-no:-1, -no:-1) = real( Le_ii(np-no:po, np-no:po,  0), RDP )
        Ld(-no:-1,   0:po) = real( Le_ii(np-no:po,     0:po,  1), RDP )

        ! lines from present element
        Ld(0:po, -no:-1) = real( Le_bc(0:po, np-no:po, -1), RDP )
        Ld(0:po,   0:po) = real( Le_bc(0:po,     0:po,  0), RDP )

        ! eigenvectors and eigenvalues .........................................

        allocate(Sd(-no:po,-no:po), Vd(-no:po))

        call SolveGeneralizedEigenproblem(Ld, Md, Vd, Sd)

        ! inject eigenvectors and eigenvalues
        S(-no:po,-no:po) = real( Sd, kind = RWP )
        V(-no:po)        = real( Vd, kind = RWP )

        ! weights ..............................................................

        W(-no:po) = real( Ws(-no:po), kind = RWP )

        ! right boundary zone
        do k = po-no+1, po
          W(k) = W(k) + real(Ws(2*po + 1 - k), RWP)
        end do

      end if

    end associate

  end subroutine ComputeGenericOperators

  !-----------------------------------------------------------------------------
  !> Set subdomain configurations and eigenvalues: const isotropic w/o SVV

  subroutine InitSubdomains_CI(this, mesh, lambda, nu, bc)
    class(DG_SchwarzOperator_3D_RWP), intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh   !< mesh partition
    real(RNP),      intent(in) :: lambda !< Helmholtz parameter
    real(RNP),      intent(in) :: nu     !< physical diffusivity
    character,      intent(in) :: bc(:)  !< BC {'D','N','P'}

    call InitSubdomains_CS(this, mesh, lambda, nu, ZERO, bc)

  end subroutine InitSubdomains_CI

  !-----------------------------------------------------------------------------
  !> Set subdomain configurations and eigenvalues: const isotropic with SVV

  subroutine InitSubdomains_CS(this, mesh, lambda, nu, nu_svv, bc)
    class(DG_SchwarzOperator_3D_RWP), intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh   !< mesh partition
    real(RNP),      intent(in) :: lambda !< Helmholtz parameter
    real(RNP),      intent(in) :: nu     !< physical diffusivity
    real(RNP),      intent(in) :: nu_svv !< spectral diffusivity
    character,      intent(in) :: bc(:)  !< BC {'D','N','P'}

    character :: bc_face(size(bc))
    real(RNP) :: dx(3), dy(3), dz(3)
    real(RNP) :: l1, l2, l3, lmb
    real(RNP) :: g0, g1, g2, g3, nu_0
    integer   :: e, i, j, k, ns
    integer   :: cfg(3)

    ! preliminaries ............................................................

    ns = size(this % V, 1)

    !$omp master

    if (allocated(this % cfg  )) deallocate(this % cfg  )
    if (allocated(this % D_inv)) deallocate(this % D_inv)

    allocate(this % cfg(3, mesh%n_elem))
    allocate(this % D_inv(ns, ns, ns, mesh%n_elem))

    !$omp end master
    !$omp barrier

    nu_0 = nu + nu_svv

    ! cfg and D_inv ............................................................

    associate(V => this % V, W => this % W, D_inv => this % D_inv)

      !$omp do
      do e = 1, mesh % n_elem

        ! get element face boundary conditions
        do i = 1, 6
          j = mesh % element(e) % face(i) % boundary
          if (j > 0) then
            bc_face(i) = bc(j)
          else
            bc_face(i) = ''
          end if
        end do

        ! identify subdomain configuration
        cfg(1) = this % ConfigurationID( bc_face(1:2) )
        cfg(2) = this % ConfigurationID( bc_face(3:4) )
        cfg(3) = this % ConfigurationID( bc_face(5:6) )
        this % cfg(:,e) = cfg

        ! extensions of the corresponding cuboid
        dx = 2 * mesh % x_cube(:,e,1) ! dx(i) = ∂x/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dy = 2 * mesh % x_cube(:,e,2) ! dy(i) = ∂y/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dz = 2 * mesh % x_cube(:,e,3) ! dz(i) = ∂z/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        l1 = sqrt(dx(1)**2 + dy(1)**2 + dz(1)**2) ! ξ₁ = ξ  extension
        l2 = sqrt(dx(2)**2 + dy(2)**2 + dz(2)**2) ! ξ₂ = η  extension
        l3 = sqrt(dx(3)**2 + dy(3)**2 + dz(3)**2) ! ξ₃ = ζ  extension

        ! metric coefficients
        g0 = l1 * l2 * l3
        g1 = l2 * l3 / l1
        g2 = l3 * l1 / l2
        g3 = l1 * l2 / l3

        do k = 1, ns
        do j = 1, ns
        do i = 1, ns

          if (W(i,cfg(1)) > 0 .and. W(j,cfg(2)) > 0 .and. W(k,cfg(3)) > 0) then

            ! eigenvalues of approximate 3D subproblem
            lmb = lambda * g0  + nu_0 * ( g1 * V(i, cfg(1)) &
                                        + g2 * V(j, cfg(2)) &
                                        + g3 * V(k, cfg(3)) )

            D_inv(i,j,k,e) = real(ONE / lmb, kind = RWP)

          else

            D_inv(i,j,k,e) = 0

          end if

        end do
        end do
        end do

      end do

    end associate

  end subroutine InitSubdomains_CS

  !-----------------------------------------------------------------------------
  !> Set subdomain configurations and eigenvalues: variable isotropic

  subroutine InitSubdomains_VI(this, eop, mesh, lambda, nu, bc)
    class(DG_SchwarzOperator_3D_RWP), intent(inout) :: this
    class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
    class(Mesh_3D), intent(in) :: mesh               !< mesh partition
    real(RNP),      intent(in) :: lambda             !< Helmholtz parameter
    real(RNP),      intent(in) :: nu(0:,0:,0:,:)     !< physical diffusivity
    character,      intent(in) :: bc(:)              !< BC {'D','N','P'}

    character :: bc_face(size(bc))
    real(RNP) :: A(0:eop%po, 0:eop%po, 0:eop%po)
    real(RNP) :: dx(3), dy(3), dz(3)
    real(RNP) :: l1, l2, l3, lmb
    real(RNP) :: g0, g1, g2, g3, nu_0
    integer   :: e, i, j, k, ns, po
    integer   :: cfg(3)

    ! preliminaries ............................................................

    ns = size(this % V, 1)
    po = eop % po

    !$omp master

    if (allocated(this % cfg  )) deallocate(this % cfg  )
    if (allocated(this % D_inv)) deallocate(this % D_inv)

    allocate(this % cfg(3, mesh%n_elem))
    allocate(this % D_inv(ns, ns, ns, mesh%n_elem))

    !$omp end master
    !$omp barrier

    ! averaging operator
    do k = 0, po
    do j = 0, po
    do i = 0, po
      A(i,j,k) = ONE/8 * eop%w(i) * eop%w(j) * eop%w(k)
    end do
    end do
    end do

    ! cfg and D_inv ............................................................

    associate(V => this % V, W => this % W, D_inv => this % D_inv)

      !$omp do
      do e = 1, mesh % n_elem

        ! get element face boundary conditions
        do i = 1, 6
          j = mesh % element(e) % face(i) % boundary
          if (j > 0) then
            bc_face(i) = bc(j)
          else
            bc_face(i) = ''
          end if
        end do

        ! identify subdomain configuration
        cfg(1) = this % ConfigurationID( bc_face(1:2) )
        cfg(2) = this % ConfigurationID( bc_face(3:4) )
        cfg(3) = this % ConfigurationID( bc_face(5:6) )
        this % cfg(:,e) = cfg

        ! extensions of the corresponding cuboid
        dx = 2 * mesh % x_cube(:,e,1) ! dx(i) = ∂x/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dy = 2 * mesh % x_cube(:,e,2) ! dy(i) = ∂y/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dz = 2 * mesh % x_cube(:,e,3) ! dz(i) = ∂z/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        l1 = sqrt(dx(1)**2 + dy(1)**2 + dz(1)**2) ! ξ₁ = ξ  extension
        l2 = sqrt(dx(2)**2 + dy(2)**2 + dz(2)**2) ! ξ₂ = η  extension
        l3 = sqrt(dx(3)**2 + dy(3)**2 + dz(3)**2) ! ξ₃ = ζ  extension

        ! metric coefficients
        g0 = l1 * l2 * l3
        g1 = l2 * l3 / l1
        g2 = l3 * l1 / l2
        g3 = l1 * l2 / l3

        ! average diffusivity
        nu_0 = 0
        do k = 0, po
        do j = 0, po
        do i = 0, po
          nu_0 = nu_0 + A(i,j,k) * nu(i,j,k,e)
        end do
        end do
        end do

        do k = 1, ns
        do j = 1, ns
        do i = 1, ns

          if (W(i,cfg(1)) > 0 .and. W(j,cfg(2)) > 0 .and. W(k,cfg(3)) > 0) then

            ! eigenvalues of approximate 3D subproblem
            lmb = lambda * g0  + nu_0 * ( g1 * V(i, cfg(1)) &
                                        + g2 * V(j, cfg(2)) &
                                        + g3 * V(k, cfg(3)) )

            D_inv(i,j,k,e) = real(ONE / lmb, kind = RWP)

          else

            D_inv(i,j,k,e) = 0

          end if

        end do
        end do
        end do

      end do

    end associate

  end subroutine InitSubdomains_VI

  !=============================================================================
