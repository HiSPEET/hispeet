!> summary:  DG Schwarz operator for elliptic equations
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI
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
    integer   :: wp        = RDP   !< working precision {RDP,RSP}
    real(RNP) :: delta     = 0.125 !< relative overlap ≤ 1
    integer   :: no_min    = 1     !< min overlap in points
    integer   :: weighting = 5     !< weighting method {0,1,3,5,7,9}
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

    type(SuboperatorsDP) :: ops_dp   !< double precision operators
    type(SuboperatorsSP) :: ops_sp   !< single precision operators

  contains

    procedure         :: BuildSubdomainMetrics
    procedure, nopass :: GetSubdomainConfigurations

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
    call XMPI_Bcast(this % delta     , root, comm)
    call XMPI_Bcast(this % no_min    , root, comm)
    call XMPI_Bcast(this % weighting , root, comm)

  end subroutine Bcast_DG_SchwarzOptions_3D

  !=============================================================================
  ! Type-bound procedures of DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Constructor

  function New_DG_SchwarzOperator_3D(opt, eop, mesh) result(this)
    class(DG_SchwarzOptions_3D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    class(Mesh_3D),                intent(in) :: mesh   !< mesh partition

    type(DG_SchwarzOperator_3D) :: this

    call InitSchwarzOperator(this, opt, eop, mesh)

  end function New_DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Build the 1D eigenvalues, eigenvectors and weights

  subroutine InitSchwarzOperator(this, opt, eop, mesh)
    class(DG_SchwarzOperator_3D),  intent(inout) :: this
    class(DG_SchwarzOptions_3D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    class(Mesh_3D),                intent(in) :: mesh   !< mesh partition

    ! local variables ..........................................................

    real(RNP), allocatable :: Ws(:)
    real(RNP), allocatable :: S(:,:), V(:), W(:)
    character :: bc(2)

    integer :: nb, nc, no, ns, po
    integer :: i, j, k

    ! initialization ...........................................................

    po = eop % po                  ! polynomial order of elements
    nb = size(DG_SCHWARZ_BC_3D)    ! number of supported boundary conditions
    nc = nb ** 2                   ! number of 1D boundary configurations

    ! number of overlapped points
    no = count(eop % x <= 2 * opt%delta - 1)  ! apply overlap
    no = max(no, opt % no_min)                ! apply minimum
    no = max(0, min(no, po+1))                ! enforce bounds

    ! number of subdomain points per direction
    ns = po + 1 + 2*no

    allocate(Ws(ns), S(ns,ns), V(ns), W(ns))

    !$omp master

    this % po = po
    this % no = no
    this % nc = nc
    this % restrictive = opt % weighting == 9 .or. no == 0

    if (allocated( this % ops_dp % S )) deallocate( this % ops_dp % S )
    if (allocated( this % ops_dp % V )) deallocate( this % ops_dp % V )
    if (allocated( this % ops_dp % W )) deallocate( this % ops_dp % W )
    if (allocated( this % ops_dp % g )) deallocate( this % ops_dp % g )

    if (allocated( this % ops_sp % S )) deallocate( this % ops_sp % S )
    if (allocated( this % ops_sp % V )) deallocate( this % ops_sp % V )
    if (allocated( this % ops_sp % W )) deallocate( this % ops_sp % W )
    if (allocated( this % ops_sp % g )) deallocate( this % ops_sp % g )

    if (opt % wp == RSP) then
      this % wp = RSP
      allocate( this % ops_sp % S(ns,ns,nc), source = 0E0 )
      allocate( this % ops_sp % V(ns,nc)   , source = 1E0 )
      allocate( this % ops_sp % W(ns,nc)   , source = 0E0 )
    else
      this % wp = RDP
      allocate( this % ops_dp % S(ns,ns,nc), source = 0D0 )
      allocate( this % ops_dp % V(ns,nc)   , source = 1D0 )
      allocate( this % ops_dp % W(ns,nc)   , source = 0D0 )
    end if

    !$omp end master
    !$omp barrier

    ! eigensystems and weights .................................................

    call WeightDistribution(eop%x, this%no, opt%weighting, Ws)

    !$omp do collapse(2)
    do j = 1, nb
    do i = 1, nb

      k = i + nb * (j - 1)

      bc(1) = DG_SCHWARZ_BC_3D(i)
      bc(2) = DG_SCHWARZ_BC_3D(j)

      call eop % Get_SchwarzSuboperators(this%no, bc, Ws, S, V, W)

      if (this % wp == RDP) then
        this % ops_dp % S(:,:,k) = real(S, RDP)
        this % ops_dp % V(:,k)   = real(V, RDP)
        this % ops_dp % W(:,k)   = real(W, RDP)
      else
        this % ops_sp % S(:,:,k) = real(S, RSP)
        this % ops_sp % V(:,k)   = real(V, RSP)
        this % ops_sp % W(:,k)   = real(W, RSP)
      end if

    end do
    end do

    call this % BuildSubdomainMetrics(mesh)

  end subroutine InitSchwarzOperator

  !-----------------------------------------------------------------------------
  !> Build metric coefficients of subdomains

  subroutine BuildSubdomainMetrics(this, mesh)
    class(DG_SchwarzOperator_3D), intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh   !< mesh partition

    real(RNP) :: dx(3)
    integer   :: e

    !$omp master
    if (allocated(this % ops_sp % g)) deallocate(this % ops_sp % g)
    if (allocated(this % ops_dp % g)) deallocate(this % ops_dp % g)
    if (this % wp == RSP) then
      allocate(this % ops_sp % g(4, mesh%n_elem_active))
    else
      allocate(this % ops_dp % g(4, mesh%n_elem_active))
    end if
    !$omp end master
    !$omp barrier

    !$omp do
    do e = 1, mesh % n_elem_active

      ! extensions of the corresponding cuboid
      call mesh % element(e) % GetCuboidDimensions(dx)

      ! metric coefficients
      if (this % wp == RSP) then
        this % ops_sp % g(1,e) = real(dx(2) * dx(3) / dx(1) , RSP)
        this % ops_sp % g(2,e) = real(dx(3) * dx(1) / dx(2) , RSP)
        this % ops_sp % g(3,e) = real(dx(1) * dx(2) / dx(3) , RSP)
        this % ops_sp % g(4,e) = real(dx(1) * dx(2) * dx(3) , RSP)
      else
        this % ops_dp % g(1,e) = real(dx(2) * dx(3) / dx(1) , RDP)
        this % ops_dp % g(2,e) = real(dx(3) * dx(1) / dx(2) , RDP)
        this % ops_dp % g(3,e) = real(dx(1) * dx(2) / dx(3) , RDP)
        this % ops_dp % g(4,e) = real(dx(1) * dx(2) * dx(3) , RDP)
      end if

    end do

  end subroutine BuildSubdomainMetrics

  !-----------------------------------------------------------------------------
  !> Get subdomain boundary configurations and metrics

  subroutine GetSubdomainConfigurations(mesh, bc, cfg)
    class(Mesh_3D), intent(in)  :: mesh     !< mesh partition
    character,      intent(in)  :: bc(:)    !< BC {'D','N','P'}
    integer,        intent(out) :: cfg(:,:) !< subdomain configurations

    character :: bc_face(6)
    integer   :: e, i, j

    !$omp do
    do e = 1, mesh % n_elem_active

      ! get element face boundary conditions
      do i = 1, 6
        j = mesh % element(e) % face(i) % boundary
        if (j > 0) then
          bc_face(i) = bc(j)
        else if (j == 0) then
          ! border to frozen element treated like Dirichlet boundary
          bc_face(i) = 'D'
        else
          bc_face(i) = ''
        end if
      end do

       cfg(1,e) = ConfigurationID( bc_face(1:2) )
       cfg(2,e) = ConfigurationID( bc_face(3:4) )
       cfg(3,e) = ConfigurationID( bc_face(5:6) )

    end do

  end subroutine GetSubdomainConfigurations

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
