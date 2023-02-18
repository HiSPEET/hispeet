!> summary:  DG Schwarz operator for elliptic equations
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use Eigenproblems, only: SolveGeneralizedEigenproblem
  use Schwarz_Weighting
  use DG__Element_Operators__1D
  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: DG_SCHWARZ_BC_3D
  public :: DG_SchwarzOptions_3D
  public :: DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Supported boundary conditions
  !>
  !> Periodic boundaries are treated as interior: 'P' → ' '

  character, parameter :: DG_SCHWARZ_BC_3D(3) = [ ' ', 'D', 'N']

  !-----------------------------------------------------------------------------
  !> Options for initializing the Schwarz operator

  type DG_SchwarzOptions_3D
    integer :: wp = RDP       !< working precision {RDP,RSP}
    integer :: no = 1         !< number of overlapped points
    integer :: weighting = 5  !< weighting method {0,1,3,5,7,9}
  contains
    procedure :: Bcast => Bcast_DG_SchwarzOptions_3D
  end type DG_SchwarzOptions_3D

  !-----------------------------------------------------------------------------
  !> Schwarz suboperators with single precision

  type SuboperatorsSP
    real(RSP), allocatable :: S(:,:,:) !< eigenvectors per config
    real(RSP), allocatable :: V(:,:)   !< eigenvalues  per config
    real(RSP), allocatable :: W(:,:)   !< weights      per config
    real(RSP), allocatable :: g(:,:)   !< metric coefficients g(4,e)
  end type SuboperatorsSP

  !-----------------------------------------------------------------------------
  !> Schwarz suboperators with double precision

  type SuboperatorsDP
    real(RDP), allocatable :: S(:,:,:) !< eigenvectors per config
    real(RDP), allocatable :: V(:,:)   !< eigenvalues  per config
    real(RDP), allocatable :: W(:,:)   !< weights      per config
    real(RDP), allocatable :: g(:,:)   !< metric coefficients g(4,e)
  end type SuboperatorsDP

  !-----------------------------------------------------------------------------
  !> Schwarz operator
  !>
  !> In the Schwarz method we consider a rectangular subdomain surrounding an
  !> element located in its center. Curvilinear elements are approximated by
  !> their cuboids. The corresponding subdomains are constructed by adopting
  !> `no` layers of collocation points from the adjoining elements.
  !>
  !> The Schwarz operator is the inverse of the truncated diffusion operator,
  !> which is given in tensor-product form by
  !>
  !>     A  =  g₁ M x M x L (ν + νˢ)
  !>        +  g₂ M x L x M (ν + νˢ)
  !>        +  g₃ L x M x M (ν + νˢ)
  !>        +  g₄ M x M x M λ
  !>
  !> where `M` is the 1D mass matrix, `L` the corresponding stiffness matrix,
  !> λ the Helmholtz parameter and ν, νˢ the physical and spectral diffusivity
  !> coefficients with a fixed ratio νˢ/(ν + νˢ).  The operators `M` and `L`
  !> are normalized to unit mesh spacing and depend on the following parameters
  !>
  !>   -  polynomial order
  !>   –  further discretization parameters such as penalty factors
  !>   -  number of overlapped points `no`
  !>   -  boundary conditions.
  !>
  !> The element extensions `∆x` are incorporated into the metric coefficients
  !>
  !>     g₁ = ∆x₂ * ∆x₃ / ∆x₁
  !>     g₂ = ∆x₃ * ∆x₁ / ∆x₂
  !>     g₃ = ∆x₁ * ∆x₂ / ∆x₃
  !>     g₄ = ∆x₁ * ∆x₂ * ∆x₃
  !>
  !> The element-boundary configuration `cfg(d,e)` describes the conditions at
  !> the element faces for each standard direction `d`:
  !>
  !>   -  In the standard configuration, the element is completely enclosed by
  !>      adjoining elements and, hence, every everywhere coated by the layer
  !>      of overlapped points.
  !>
  !>   -  In boundary configurations, one or more element faces coincide with
  !>      the boundary of the computational domain. At those faces, the
  !>      Helmholtz operator is modified according to the boundary conditions,
  !>      and no exterior points are adopted to the Schwarz subdomain.
  !>
  !> Considering interior (I), Dirichlet (D) and Neumann (N) faces, 9 different
  !> configurations have to be distinguished in each coordinate direction:
  !>
  !>    1.  I-I
  !>    2.  D-I
  !>    3.  N-I
  !>    4.  I-D
  !>    5.  D-D
  !>    6.  N-D
  !>    7.  I-N
  !>    8.  D-N
  !>    9.  N-N
  !>
  !> For computing the inverse operator, a generalized 1D eigenvalue problem is
  !> solved for every configuration in each coordinate direction, yielding the
  !> matrices of right eigenvectors `S₁`, `S₂`, `S₃` and the diagonal matrices
  !> of eigenvalues `Λ₁`, `Λ₂`, `Λ₃` such that
  !>
  !>     S₁ᵀ L₁ S₁ = Λ₁
  !>     S₁ᵀ M₁ S₁ = I₁
  !>
  !> where `I₁` ist the matching unit matrix etc. The inverse Helmholtz operator
  !> can be then expressed in the tensor-product form
  !>
  !>     A⁻¹  =  (S₃ x S₂ x S₁ D⁻¹ (S₃ᵀ x S₂ᵀ x S₁ᵀ)
  !>
  !> with the diagonal matrix
  !>
  !>     D  =  g₁ I₃ x I₂ x Λ₁ (ν + νˢ)
  !>        +  g₂ I₃ x Λ₂ x I₁ (ν + νˢ)
  !>        +  g₃ Λ₃ x I₂ x I₁ (ν + νˢ)
  !>        +  g₄ I₃ x I₂ x I₁ λ
  !>
  !> As λ, ν and νˢ can vary, `D⁻¹` is not stored, but computed on the fly using
  !> the metric coefficients and the 1D eigenvalues, which are kept within the
  !> operator.
  !>
  !> Before assembling the global correction to an approximate solution, the
  !> subdomain correction is weighted according to
  !>
  !>     Δu = W (A⁻¹ r)
  !>
  !> The weights form a diagonal tensor-product matrix
  !>
  !>     W = W₃ x W₂ x W₁
  !>
  !> with 1D distributions `W₁`, `W₂`, `W₃` depending on the element-boundary
  !> configuration.

  type DG_SchwarzOperator_3D

    integer :: wp = -1     !< working precision
    integer :: po = -1     !< polynomial order
    integer :: no = -1     !< number of overlapped layers
    integer :: nc = -1     !< number of 1D configurations
    logical :: restrictive !< T/F for ex/including neighbor results

    integer, allocatable :: cfg(:,:) !< subdomain configurations (nc,*)
    type(SuboperatorsDP) :: ops_dp   !< double precision operators
    type(SuboperatorsSP) :: ops_sp   !< single precision operators

  contains

    generic :: RestrictResidual => RestrictResidual_RDP, RestrictResidual_RSP
    procedure, private :: RestrictResidual_RDP, RestrictResidual_RSP

    generic :: MergeCorrections => MergeCorrections_RDP, MergeCorrections_RSP
    procedure, private :: MergeCorrections_RDP ,MergeCorrections_RSP

  end type DG_SchwarzOperator_3D

  ! constructor interface
  interface DG_SchwarzOperator_3D
    module procedure New_DG_SchwarzOperator_3D
  end interface

  !=============================================================================
  ! Interfaces to double precision submodule procedures

  interface

    module subroutine Restrict_Structured_RDP(this, mesh, buf_r, r, rs, sgn)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
      real(RDP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
      integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]
    end subroutine Restrict_Structured_RDP

    module subroutine Restrict_Unstructured_RDP(this, mesh, buf_r, r, rs, sgn)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
      real(RDP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
      integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]
    end subroutine Restrict_Unstructured_RDP

    module subroutine Merge_Core_RDP(this, mesh, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RDP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Core_RDP

    module subroutine Merge_Structured_RDP(this, mesh, buf_us, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RDP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Structured_RDP

    module subroutine Merge_Unstructured_RDP(this, mesh, buf_us, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RDP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Unstructured_RDP

  end interface

  !=============================================================================
  ! Interfaces to single precision submodule procedures

  interface

    module subroutine Restrict_Structured_RSP(this, mesh, buf_r, r, rs, sgn)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
      real(RSP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
      integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]
    end subroutine Restrict_Structured_RSP

    module subroutine Restrict_Unstructured_RSP(this, mesh, buf_r, r, rs, sgn)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
      real(RSP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
      integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]
    end subroutine Restrict_Unstructured_RSP

    module subroutine Merge_Core_RSP(this, mesh, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RSP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Core_RSP

    module subroutine Merge_Structured_RSP(this, mesh, buf_us, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RSP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Structured_RSP

    module subroutine Merge_Unstructured_RSP(this, mesh, buf_us, us, u)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
      real(RSP),      intent(inout) :: us(:,:,:,:)   !< mesh variable
    end subroutine Merge_Unstructured_RSP

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of DG_SchwarzOptions_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type DG_SchwarzOptions_3D

   subroutine Bcast_DG_SchwarzOptions_3D(this, root, comm)
    class(DG_SchwarzOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % wp        , root, comm)
    call XMPI_Bcast(this % no        , root, comm)
    call XMPI_Bcast(this % weighting , root, comm)

  end subroutine Bcast_DG_SchwarzOptions_3D

  !=============================================================================
  ! Type-bound procedures of DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Constructor

  function New_DG_SchwarzOperator_3D(opt, eop, mesh, bc, r_nu_s) result(this)
    class(DG_SchwarzOptions_3D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    class(Mesh_3D),                intent(in) :: mesh   !< mesh partition
    character,                     intent(in) :: bc(:)  !< BC {'D','N','P'}
    real(RNP),           optional, intent(in) :: r_nu_s !< ratio νˢ/(ν + νˢ) [0]

    type(DG_SchwarzOperator_3D) :: this

    call InitSchwarzOperator (this, opt, eop, mesh, r_nu_s)
    call InitSchwarzDomains  (this, mesh, bc)

  end function New_DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Build the 1D eigenvalues, eigenvectors and weights

  subroutine InitSchwarzOperator(this, opt, eop, mesh, r_nu_s)
    class(DG_SchwarzOperator_3D),  intent(inout) :: this
    class(DG_SchwarzOptions_3D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    class(Mesh_3D),                intent(in) :: mesh   !< mesh partition
    real(RNP),           optional, intent(in) :: r_nu_s !< ratio νˢ/(ν+νˢ) [0]

    ! local variables ..........................................................

    real(RNP), allocatable :: Ws(:)
    real(RDP), allocatable :: S(:,:), V(:), W(:)
    character :: bc(2)

    integer :: nb, nc, ne, no, ns, po
    integer :: i, j, k

    ! initialization ...........................................................

    po = eop % po                  ! polynomial order of elements
    nb = size(DG_SCHWARZ_BC_3D)    ! number of supported boundary conditions
    nc = nb ** 2                   ! number of 1D boundary configurations
    ne = mesh % n_elem             ! number of local elements / subdomains
    no = max(0, min(opt%no, po+1)) ! number of overlapped points
    ns = po + 1 + 2*no             ! number of subdomain points per direction

    allocate(Ws(ns), S(ns,ns), V(ns), W(ns))

    !$omp single

    this % po = po
    this % no = no
    this % nc = nc
    this % restrictive = opt % weighting == 9 .or. no == 0

    if (allocated(this % cfg)) deallocate(this % cfg)

    if (allocated( this % ops_dp % S )) deallocate( this % ops_dp % S )
    if (allocated( this % ops_dp % V )) deallocate( this % ops_dp % V )
    if (allocated( this % ops_dp % W )) deallocate( this % ops_dp % W )
    if (allocated( this % ops_dp % g )) deallocate( this % ops_dp % g )

    if (allocated( this % ops_sp % S )) deallocate( this % ops_sp % S )
    if (allocated( this % ops_sp % V )) deallocate( this % ops_sp % V )
    if (allocated( this % ops_sp % W )) deallocate( this % ops_sp % W )
    if (allocated( this % ops_sp % g )) deallocate( this % ops_sp % g )

    allocate(this % cfg(3,ne))

    if (opt % wp == RSP) then
      this % wp = RSP
      allocate( this % ops_sp % S(ns,ns,nc), source = 0E0 )
      allocate( this % ops_sp % V(ns,nc)   , source = 1E0 )
      allocate( this % ops_sp % W(ns,nc)   , source = 0E0 )
      allocate( this % ops_sp % g(4,ne)                   )
    else
      this % wp = RDP
      allocate( this % ops_dp % S(ns,ns,nc), source = 0D0 )
      allocate( this % ops_dp % V(ns,nc)   , source = 1D0 )
      allocate( this % ops_dp % W(ns,nc)   , source = 0D0 )
      allocate( this % ops_dp % g(4,ne)                   )
    end if

    !$omp end single

    ! eigensystems and weights .................................................

    call WeightDistribution(eop%x, this%no, opt%weighting, Ws)

    !$omp do collapse(2)
    do j = 1, nb
    do i = 1, nb

      k = i + nb * (j - 1)

      bc(1) = DG_SCHWARZ_BC_3D(i)
      bc(2) = DG_SCHWARZ_BC_3D(j)

      call InitSuboperators(eop, this%no, bc, Ws, S, V, W, r_nu_s)

      if (this % wp == RDP) then
        this % ops_dp % S(:,:,k) = S
        this % ops_dp % V(:,k)   = V
        this % ops_dp % W(:,k)   = W
      else
        this % ops_sp % S(:,:,k) = real(S, RSP)
        this % ops_sp % V(:,k)   = real(V, RSP)
        this % ops_sp % W(:,k)   = real(W, RSP)
      end if

    end do
    end do

  end subroutine InitSchwarzOperator

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
  !> entries are retained for regularity. They are, however set to `0` in the
  !> eigenvectors and weights. The corresponding eigenvalues are set `1` in
  !> order to avoid floating points exceptions when used as a divisor.


  subroutine InitSuboperators(eop, no, bc, Ws, S, V, W, r_nu_s)
    class(DG_ElementOperators_1D), intent(in) :: eop !< IP-DG element operators
    integer,   intent(in)  :: no              !< overlap
    character, intent(in)  :: bc(2)           !< left/right boundary conditions
    real(RNP), intent(in)  :: Ws(-no:)        !< standard weights
    real(RDP), intent(out) :: S(-no:,-no:)    !< subdomain eigenvectors
    real(RDP), intent(out) :: V(-no:)         !< subdomain eigenvalues
    real(RDP), intent(out) :: W(-no:)         !< subdomain weights
    real(RNP), optional, intent(in) :: r_nu_s !< ratio νˢ/(ν + νˢ) [0]

    character, parameter   :: ii(2) = [ ' ', ' ' ]
    real(RNP), parameter   :: dx(-1:1) = 1
    real(RNP), allocatable :: Le_ii(:,:,:), Le_bc(:,:,:)
    real(RDP), allocatable :: Ld(:,:), Md(:), Sd(:,:), vd(:)

    real(RNP) :: nu, nu_s
    integer   :: po, np
    integer   :: k

    po = eop%po
    np = po + 1

    allocate(Le_ii(0:po,0:po,-1:1))
    allocate(Le_bc(0:po,0:po,-1:1))

    if (present(r_nu_s)) then
      nu_s = r_nu_s
    else
      nu_s = ZERO
    end if
    nu = max(ONE - nu_s, ZERO)

    ! stiffness matrix for interior element with dx=1, bc = ii
    call eop % Get_DiffusionMatrix(dx, ii, nu, nu_s, Ae = Le_ii)

    ! stiffness matrix for given boundary conditions
    call eop % Get_DiffusionMatrix(dx, bc, nu, nu_s, Ae = Le_bc)

    ! initialization of eigenvectors S, eigenvalues V and weights W
    S = 0
    V = 1
    W = 0

    associate(Ms => eop % w)

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
        S = Sd
        V = Vd

        ! weights ..............................................................

        W = real( Ws, kind = RDP )

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
        S(0:po,0:po) = Sd
        V(0:po)      = Vd

        ! weights ..............................................................

        W(0:po) = real( Ws(0:po), kind = RDP )

        ! left boundary zone
        do k = 0, no-1
          W(k) = W(k) + real(Ws(-k-1), RDP)
        end do

        ! right boundary zone
        do k = po-no+1, po
          W(k) = W(k) + real(Ws(2*po + 1 - k), RDP)
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

        S(0:po+no,0:po+no) = Sd
        V(0:po+no)         = Vd

        ! weights ..............................................................

        W(0:) = real( Ws(0:), kind = RDP )

        ! boundary zone
        do k = 0, no-1
          W(k) = W(k) + real(Ws(-k-1), RDP)
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
        S(-no:po,-no:po) = Sd
        V(-no:po)        = Vd

        ! weights ..............................................................

        W(-no:po) = real( Ws(-no:po), kind = RDP )

        ! right boundary zone
        do k = po-no+1, po
          W(k) = W(k) + real(Ws(2*po + 1 - k), RDP)
        end do

      end if

    end associate

  end subroutine InitSuboperators

  !-----------------------------------------------------------------------------
  !> Build subdomain configurations and metrics

  subroutine InitSchwarzDomains(this, mesh, bc)
    class(DG_SchwarzOperator_3D), intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh   !< mesh partition
    character,      intent(in) :: bc(:)  !< BC {'D','N','P'}

    character :: bc_face(6)
    real(RNP) :: dx(3), dy(3), dz(3)
    real(RNP) :: l1, l2, l3
    integer   :: c1, c2, c3
    integer   :: e, i, j

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

      ! set subdomain configuration
      c1 = ConfigurationID( bc_face(1:2) )
      c2 = ConfigurationID( bc_face(3:4) )
      c3 = ConfigurationID( bc_face(5:6) )
      this % cfg(1,e) = c1
      this % cfg(2,e) = c2
      this % cfg(3,e) = c3

      ! extensions of the corresponding cuboid
      associate(x_c => mesh % element(e) % geometry % x_c)
        dx = 2 * x_c(1:3,1) ! dx(i) = ∂x/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dy = 2 * x_c(1:3,2) ! dy(i) = ∂y/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        dz = 2 * x_c(1:3,3) ! dz(i) = ∂z/∂ξᵢ ∆ξᵢ  with  ∆ξᵢ = 2
        l1 = sqrt(dx(1)**2 + dy(1)**2 + dz(1)**2) ! ξ₁ = ξ  extension
        l2 = sqrt(dx(2)**2 + dy(2)**2 + dz(2)**2) ! ξ₂ = η  extension
        l3 = sqrt(dx(3)**2 + dy(3)**2 + dz(3)**2) ! ξ₃ = ζ  extension
      end associate

      ! metric coefficients
      if (this % wp == RSP) then
        this % ops_sp % g(1,e) = real(l2 * l3 / l1 , RSP)
        this % ops_sp % g(2,e) = real(l3 * l1 / l2 , RSP)
        this % ops_sp % g(3,e) = real(l1 * l2 / l3 , RSP)
        this % ops_sp % g(4,e) = real(l1 * l2 * l3 , RSP)
      else
        this % ops_dp % g(1,e) = real(l2 * l3 / l1 , RDP)
        this % ops_dp % g(2,e) = real(l3 * l1 / l2 , RDP)
        this % ops_dp % g(3,e) = real(l1 * l2 / l3 , RDP)
        this % ops_dp % g(4,e) = real(l1 * l2 * l3 , RDP)
      end if

    end do

  end subroutine InitSchwarzDomains

  !-----------------------------------------------------------------------------
  !> Restrict mesh variable to subdomains -- double precision
  !>
  !> The mesh variable must be dimensioned as `r(0:po,0:po,0:po,ne+ng)`, where
  !> `ne` is the number of local elements and `ng` the number of ghosts.

  subroutine RestrictResidual_RDP(this, mesh, buf_r, r, rs, sgn)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
    real(RDP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
    integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]

    if (mesh % structured) then
      call Restrict_Structured_RDP(this, mesh, buf_r, r, rs, sgn)
    else
      call Restrict_Unstructured_RDP(this, mesh, buf_r, r, rs, sgn)
    end if

  end subroutine RestrictResidual_RDP

  !-----------------------------------------------------------------------------
  !> Merge subdomain contributions into mesh variable -- double precision
  !>
  !> The subdomain variable must be dimensioned as `u(ns,ns,ns,ne+ng)`, where
  !> `ns` is the numbeer of subdomain points per direction, `ne` the number of
  !> local elements and `ng` the number of ghosts.

  subroutine MergeCorrections_RDP(this, mesh, buf_us, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
    real(RDP),      intent(inout) :: us(:,:,:,:)   !< mesh variable

    if (this % restrictive) then ! merge core regions only
      call Merge_Core_RDP(this, mesh, us, u)
    else if (mesh % structured) then
      call Merge_Structured_RDP(this, mesh, buf_us, us, u)
    else
      call Merge_Unstructured_RDP(this, mesh, buf_us, us, u)
    end if

  end subroutine MergeCorrections_RDP

  !-----------------------------------------------------------------------------
  !> Restrict mesh variable to subdomains -- single precision

  subroutine RestrictResidual_RSP(this, mesh, buf_r, r, rs, sgn)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
    real(RSP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
    integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]

    if (mesh % structured) then
      call Restrict_Structured_RSP(this, mesh, buf_r, r, rs, sgn)
    else
      call Restrict_Unstructured_RSP(this, mesh, buf_r, r, rs, sgn)
    end if

  end subroutine RestrictResidual_RSP

  !-----------------------------------------------------------------------------
  !> Merge subdomain contributions into mesh variable -- single precision

  subroutine MergeCorrections_RSP(this, mesh, buf_us, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
    real(RSP),      intent(inout) :: us(:,:,:,:)   !< mesh variable

    if (this % restrictive) then ! merge core regions only
      call Merge_Core_RSP(this, mesh, us, u)
    else if (mesh % structured) then
      call Merge_Structured_RSP(this, mesh, buf_us, us, u)
    else
      call Merge_Unstructured_RSP(this, mesh, buf_us, us, u)
    end if

  end subroutine MergeCorrections_RSP

  !=============================================================================
  ! Utilities

  !-----------------------------------------------------------------------------
  !> Returns the 1D subdomain configuration ID corresponding to the given BCs

  pure integer function ConfigurationID(bc) result(cfg)
    character, intent(in) :: bc(2) !< left/right boundary types {' ','D','N','P'}

    integer :: i1, i2

    select case(bc(1))
    case('D')
      i1 = 2
    case('N')
      i1 = 3
    case default ! ' ' and 'P'
      i1 = 1
    end select

    select case(bc(2))
    case('D')
      i2 = 2
    case('N')
      i2 = 3
    case default ! ' ' and 'P'
      i2 = 1
    end select

    cfg = i1 + (i2 - 1) * size(DG_SCHWARZ_BC_3D)

  end function ConfigurationID

  !=============================================================================

end module DG__Schwarz_Operator__3D
