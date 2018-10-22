!> summary:  Schwarz method for elliptic equation with equidistant DG
!> author:   Joerg Stiller
!> date:     2017/09/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Schwarz method for elliptic equation with equidistant DG
!>
!> @todo
!>
!>   * Explore potential run time savings by
!>
!>       - reusing element transfer buffers
!>       - reusing element-boundary configuration
!>       - reusing workspace
!>
!>   * Put restriction to and merging from subdomains in submodule(s)
!>
!> @endtodo
!===============================================================================

module CART__DG_Elliptic_CI_Schwarz

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF
  use Array_Assignments, only: AssignScalar
  use Array_Reductions,  only: ScalarProduct
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use Execution_Control, only: Error
  use Schwarz_Weighting
  use DG_Element_Operators_1D

  use XMPI

  use CART__TPO_Schwarz_Iso
  use CART__Mesh_Element
  use CART__Mesh_Partition
  use CART__Element_Transfer_Buffer
  use CART__DG_Element_Operators
  use CART__DG_Elliptic_CI_Residual

  implicit none
  private

  public :: SchwarzOperator

  !-----------------------------------------------------------------------------
  !> Schwarz operator and related procedures
  !>
  !> In the Schwarz method we consider a rectangular subdomain surrounding an
  !> element located in its center. The subdomain is constructed by adopting a
  !> layer of collocation points from the adjoining elements. The thickness of
  !> this layer is a directional property, which depends on two parameters:
  !> the relative thickness `delta` and the minimal number of overlapped points
  !> `no_min`. Typically, the thickness assumes a value between 0 and 1, though
  !> a negative value can be chosen for restricting the subdomain to the element
  !> alone.
  !>
  !> The Schwarz operator is the inverse of the truncated elliptic operator,
  !> which is given in tensor-product form by
  !>
  !>       A  =  c0 M3 x M2 x M1
  !>          +  c1 M3 x M2 x L1
  !>          +  c2 M3 x L2 x M1
  !>          +  c3 L3 x M2 x M1
  !>
  !> where `M1`, `M2`, `M3` are the 1D mass matrices and  `L1`, `L2`, `L3` the
  !> corresponding stiffness matrices. These operators are normalized zo unit
  !> mesh spacing and, thus, depend only on
  !>
  !>   *  the polynomial order and the penalty factor, both contained in the
  !>      element operators stored in `eop`,
  !>
  !>   *  the number of overlapped points `no`, and
  !>
  !>   *  the element-boundary configuration.
  !>
  !> The effect of element extensions `dx` is incorporated into the coefficients
  !>
  !>     c0 = dx(1) * dx(2) * dx(3) * lambda
  !>     c1 = dx(2) * dx(3) / dx(1) * nu
  !>     c2 = dx(3) * dx(1) / dx(2) * nu
  !>     c3 = dx(1) * dx(2) / dx(3) * nu
  !>
  !> `lambda` represents the Helmholtz parameter and  `nu` the diffusivity.
  !>
  !> The element-boundary configuration describes the conditions met at the
  !> element faces:
  !>
  !>   *  In the standard configuration, the element is completely enclosed by
  !>      adjoining elements and, hence, every everywhere coated by the layer
  !>      of overlapped points.
  !>
  !>   *  In boundary configurations, one or more element faces coincide with
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
  !> matrices of right eigenvectors `S1`, `S2`, `S3` and the diagonal matrices
  !> of eigenvalues `Λ1`, `Λ2`, `Λ3` such that
  !>
  !>     S1ᵀ L1 S1 = Λ1
  !>     S1ᵀ M1 S1 = I1
  !>
  !> where `I1` ist the matching unit matrix etc. Rearranging the eigenvalues
  !> into coefficient arrays
  !>
  !>     D1 = c0/3 I1 + c1 Λ1
  !>     D2 = c0/3 I2 + c2 Λ2
  !>     D3 = c0/3 I3 + c3 Λ3
  !>
  !> the inverse Helmholtz operator can be expressed in the tensor-product form
  !>
  !>     A⁻¹  =  (S3 x S2 x S1) D⁻¹ (S3ᵀ x S2ᵀ x S1ᵀ)
  !>
  !> with the diagonal matrix
  !>
  !>     D  =  I3 x I2 x D1  +  I3 x D2 x I1  +  D3 x I2 x I1
  !>
  !> The diagonals of `D1`, `D2`, `D3` are stored in `D1`, `D2`, and `D3`,
  !> respectively.
  !>
  !> Before assembling the global correction to an approximate solution, the
  !> subdomain correction is weighted according to
  !>
  !>     Δu = W (A⁻¹ r)
  !>
  !> The weights form a diagonal tensor-product matrix
  !>
  !>     W = W3 x W2 x W1
  !>
  !> with 1D distributions `W1`, `W2`, `W3` depending on the element-boundary
  !> configuration.
  !>
  !> It is worth noting, that eigenvectors, eigenvalues and weights coincide,
  !> if the number of overlapped points is identical in each direction.
  !> To benefit from possible optimizations, this quasi-isotropic case is
  !> indicated by setting component `isotropic` to true.

  type SchwarzOperator

    type(DG_ElementOperators3D) :: eop  !< DG element operators

    integer   :: no(3)  = -1            !< overlapped node layers
    logical   :: isotropic              !< switch to isotropic operator

    real(RNP), allocatable :: S1(:,:,:) !< eigenvectors for direction 1
    real(RNP), allocatable :: S2(:,:,:) !< eigenvectors for direction 2
    real(RNP), allocatable :: S3(:,:,:) !< eigenvectors for direction 3

    real(RNP), allocatable :: V1(:,:)   !< eigenvalues for direction 1 (Λ1)
    real(RNP), allocatable :: V2(:,:)   !< eigenvalues for direction 2 (Λ2)
    real(RNP), allocatable :: V3(:,:)   !< eigenvalues for direction 3 (Λ3)

    real(RNP), allocatable :: W1(:,:)   !< weights for direction 1
    real(RNP), allocatable :: W2(:,:)   !< weights for direction 2
    real(RNP), allocatable :: W3(:,:)   !< weights for direction 3

  contains

    procedure         :: New => New_SchwarzOperator
    procedure, nopass :: IdentifyConfiguration
    procedure         :: Iteration => Iteration_C

    final :: Delete_SchwarzOperator

  end type SchwarzOperator

  !=============================================================================
  ! Private module variables

  !-----------------------------------------------------------------------------
  !> List of supported one-dimensional element-boundary configurations

  character, parameter ::         &
    boundary_configuration(2,9)   &
        = reshape( [ ' ', ' ',    &! 1: Interior  - Interior
                     'D', ' ',    &! 2: Dirichlet - Interior
                     ' ', 'D',    &! 3: Interior  - Dirichlet
                     'N', ' ',    &! 4: Neumann   - Interior
                     ' ', 'N',    &! 5: Interior  - Neumann
                     'D', 'D',    &! 6: Dirichlet - Dirichlet
                     'D', 'N',    &! 7: Dirichlet - Neumann
                     'N', 'D',    &! 8: Neumann   - Dirichlet
                     'N', 'N' ],  &! 9: Neumann   - Neumann
                   [2,9] )

contains

!-------------------------------------------------------------------------------
!> Initialization of a SchwarzOperator object

subroutine New_SchwarzOperator(this, eop, delta, no_min, weighting)

  ! arguments ..................................................................

  class(SchwarzOperator),       intent(inout) :: this !< Schwarz operator
  class(DG_ElementOperators3D), intent(in)    :: eop  !< DG element operators

  real(RNP),         intent(in) :: delta(3)  !< relative overlap
  integer, optional, intent(in) :: no_min    !< min overlap in points
  integer, optional, intent(in) :: weighting !< weighting method [5]

  ! local variables ............................................................

  real(RNP), allocatable :: W1(:), W2(:), W3(:)

  integer :: n1, n2, n3, nc
  integer :: i, k

  ! clean-up ...................................................................

  call Delete_SchwarzOperator(this)

  ! adoptions ..................................................................

  this % eop = eop

  associate(po => eop%po, Ms => eop%w)

    ! overlap ..................................................................

    do i = 1, 3
      this%no(i) = count(eop%x <= 2*delta(i) - 1)
    end do

    if (present(no_min)) then
      this%no = max(this%no, min(no_min, po+1))
    end if

    this%isotropic = all(this%no(2:3) == this%no(1))

    ! preparations .............................................................

    ! dimensions
    n1 = po + 1 + 2 * this%no(1)
    n2 = po + 1 + 2 * this%no(2)
    n3 = po + 1 + 2 * this%no(3)
    nc = size(boundary_configuration, 2)

    ! weights for standard (interior-interior) configuration
    allocate(W1(n1), W2(n3), W3(n3))
    if (present(weighting)) then
      k = weighting
    else
      k = 5
    end if
    call WeightDistribution(eop%x, this%no(1), k, W1)
    call WeightDistribution(eop%x, this%no(2), k, W2)
    call WeightDistribution(eop%x, this%no(3), k, W3)

    allocate(this%S1(n1,n1,nc), this%V1(n1,nc), this%W1(n1,nc))
    allocate(this%S2(n2,n2,nc), this%V2(n2,nc), this%W2(n2,nc))
    allocate(this%S3(n3,n3,nc), this%V3(n3,nc), this%W3(n3,nc))

    ! eigenvectors, eigenvalues and weights ....................................

    do i = 1, nc

      call GetSubdomainOperators( eop,                              &
                                  no = this%no(1),                  &
                                  bc = boundary_configuration(:,i), &
                                  Ws = W1,                          &
                                  S  = this%S1(:,:,i),              &
                                  V  = this%V1(:,i),                &
                                  W  = this%W1(:,i)                 )

      if (this%isotropic) then

        this % S2 = this % S1
        this % V2 = this % V1
        this % W2 = this % W1

        this % S3 = this % S1
        this % V3 = this % V1
        this % W3 = this % W1

      else

        call GetSubdomainOperators( eop,                              &
                                    no = this%no(2),                  &
                                    bc = boundary_configuration(:,i), &
                                    Ws = W2,                          &
                                    S  = this%S2(:,:,i),              &
                                    V  = this%V2(:,i),                &
                                    W  = this%W2(:,i)                 )

        call GetSubdomainOperators( eop,                              &
                                    no = this%no(3),                  &
                                    bc = boundary_configuration(:,i), &
                                    Ws = W3,                          &
                                    S  = this%S3(:,:,i),              &
                                    V  = this%V3(:,i),                &
                                    W  = this%W3(:,i)                 )
      end if

    end do

  end associate

end subroutine New_SchwarzOperator

!-------------------------------------------------------------------------------
!> Computes the eigenvectors, eigenvalues and weights for a 1D subdomain
!>
!> This procedure provides the eigenvectors and eigenvalues for 1D subdomains
!> assuming a constant element width of dx=1. The actual element width is
!> considered be adjusting the coefficients c0, c1, c2, and c3, as detailed in
!> the description of `SchwarzOperator`.
!> The argument `bc` specifies the left and right boundary conditions.
!> `D`indicates a Dirichlet boundary and `N` a Neumann boundary.
!> Otherwise, the existence of a neighbor element is assumed.
!>
!> Though the exterior node layers are cut off at boundaries, the corresponding
!> entries are retained for regularity. They are, however, set to `0` in the
!> eigenvectors and weights. The corresponding eigenvalues they are set `1` in
!> order to avoid floating points exceptions when used as a divisor.

subroutine GetSubdomainOperators(eop, no, bc, Ws, S, V, W)
  type(DG_ElementOperators3D), intent(in) :: eop !< DG element operators
  integer,    intent(in)  :: no           !< overlap
  character,  intent(in)  :: bc(2)        !< left/right boundary conditions
  real(RNP),  intent(in)  :: Ws(-no:)     !< standard weights
  real(RNP),  intent(out) :: S(-no:,-no:) !< subdomain eigenvectors
  real(RNP),  intent(out) :: V(-no:)      !< subdomain eigenvalues
  real(RNP),  intent(out) :: W(-no:)      !< subdomain weights

  real(RNP), parameter   :: dx(-1:1) = ONE
  real(RNP), allocatable :: Le_ii(:,:,:), Le_bc(:,:,:)
  real(RNP), allocatable :: Ld(:,:), Md(:), Sd(:,:), vd(:)

  character, parameter :: ii(2) = [ ' ', ' ' ]
  integer :: po, np
  integer :: k

  po = eop%po
  np = po + 1

  allocate(Le_ii(0:po,0:po,-1:1))
  allocate(Le_bc(0:po,0:po,-1:1))

  ! stiffness matrix in direction 1 for interior element
  call eop % Get_1D_StiffnessMatrix(direction=1, bc=ii, Le=Le_ii)

  ! stiffness matrix in direction 1 for given boundary conditions
  call eop % Get_1D_StiffnessMatrix(direction=1, bc=bc, Le=Le_bc)

  ! scale stiffness matrices to spacing of dx = 1
  Le_ii = eop%dx(1) * Le_ii
  Le_bc = eop%dx(1) * Le_bc

  ! initialization: eigenvalues V set to 1 to avoid division by zero in
  ! in case of truncated overlap zones
  S = 0
  V = 1
  W = 0

  associate(Ms => eop % W)

    if (all(bc == ' ')) then

      ! interior-interior configuration ----------------------------------------

      !!! all elements regarded as interior ones !!!

      ! mass matrix ............................................................

      allocate(Md(-no:po+no))

      Md(-no:-1   ) = HALF * Ms( np-no:po   )
      Md(  0:po   ) = HALF * Ms(     0:po   )
      Md( np:po+no) = HALF * Ms(     0:no-1 )

      ! stiffness matrix .......................................................

      allocate(Ld(-no:po+no,-no:po+no), source = ZERO)

      ! lines from preceding element
      Ld(-no:-1, -no:-1) = Le_ii(np-no:po, np-no:po,  0)
      Ld(-no:-1,   0:po) = Le_ii(np-no:po,     0:po,  1)

      ! lines from present element
      Ld(0:po, -no:-1   ) = Le_ii(0:po, np-no:po  , -1)
      Ld(0:po,   0:po   ) = Le_ii(0:po,     0:po  ,  0)
      Ld(0:po,  np:po+no) = Le_ii(0:po,     0:no-1,  1)

      ! lines from following element
      Ld(np:po+no,  0:po   ) = Le_ii(0:no-1, 0:po  , -1)
      Ld(np:po+no, np:po+no) = Le_ii(0:no-1, 0:no-1,  0)

      ! eigenvectors and eigenvalues ...........................................

      allocate(Sd, mold = Ld)
      allocate(vd, mold = Md)

      call SolveGeneralizedEigenproblem(Ld, Md, vd, Sd)

      ! inject eigenvectors and eigenvalues
      S = Sd
      V = vd

      ! weights ................................................................

      W = Ws

    else if (all(bc /= ' ')) then

      ! boundary-boundary configuration-----------------------------------------

      ! mass matrix ............................................................

      allocate(Md(0:po), source = HALF * Ms)

      ! stiffness matrix .......................................................

      allocate(Ld(0:po,0:po), source = Le_bc(:,:,0))

      ! eigenvectors and eigenvalues ...........................................

      allocate(Sd(0:po,0:po), vd(0:po))

      call SolveGeneralizedEigenproblem(Ld, Md, vd, Sd)

      ! inject eigenvectors and eigenvalues
      S(0:po,0:po) = Sd
      V(0:po)      = vd

      ! weights ................................................................

      W(0:po) = Ws(0:po)

      ! left boundary zone
      do k = 0, no-1
        W(k) = W(k) + Ws(-k-1)
      end do

      ! right boundary zone
      do k = po-no+1, po
        W(k) = W(k) + Ws(2*po + 1 - k)
      end do

    else if (bc(1) /= ' ') then

      ! boundary-interior configuration-----------------------------------------

      !!! neighbor on the right always regarded as an interior element !!!

      ! mass matrix ............................................................

      allocate(Md(0:po+no))

      Md(  0:po   ) = HALF * Ms(     0:po   )
      Md( np:po+no) = HALF * Ms(     0:no-1 )

      ! stiffness matrix .......................................................

      allocate(Ld(0:po+no,0:po+no), source = ZERO)

      ! lines from present element
      Ld( 0:po,  0:po   ) = Le_bc(0:po, 0:po  ,  0)
      Ld( 0:po, np:po+no) = Le_bc(0:po, 0:no-1,  1)

      ! lines from following element
      Ld(np:po+no,  0:po   ) = Le_ii(0:no-1, 0:po  , -1)
      Ld(np:po+no, np:po+no) = Le_ii(0:no-1, 0:no-1,  0)

      ! eigenvectors and eigenvalues ...........................................

      allocate(Sd(0:po+no,0:po+no), vd(0:po+no))

      call SolveGeneralizedEigenproblem(Ld, Md, vd, Sd)

      ! inject eigenvectors and eigenvalues

      S(0:po+no,0:po+no) = Sd
      V(0:po+no)         = vd

      ! weights ................................................................

      W(0:) = Ws(0:)

      ! boundary zone
      do k = 0, no-1
        W(k) = W(k) + Ws(-k-1)
      end do

    else

      ! interior-boundary configuration-----------------------------------------

      !!! neighbor on the left always regarded as an interior element !!!

      ! mass matrix ............................................................

      allocate(Md(-no:po))

      Md(-no:-1   ) = HALF * Ms( np-no:po   )
      Md(  0:po   ) = HALF * Ms(     0:po   )

      ! stiffness matrix .......................................................

      allocate(Ld(-no:po,-no:po), source = ZERO)

      ! lines from preceding element
      Ld(-no:-1, -no:-1) = Le_ii(np-no:po, np-no:po,  0)
      Ld(-no:-1,   0:po) = Le_ii(np-no:po,     0:po,  1)

      ! lines from present element
      Ld(0:po, -no:-1) = Le_bc(0:po, np-no:po, -1)
      Ld(0:po,   0:po) = Le_bc(0:po,     0:po,  0)

      ! eigenvectors and eigenvalues ...........................................

      allocate(Sd(-no:po,-no:po), vd(-no:po))

      call SolveGeneralizedEigenproblem(Ld, Md, vd, Sd)

      ! inject eigenvectors and eigenvalues
      S(-no:po,-no:po) = Sd
      V(-no:po)        = vd

      ! weights ................................................................

      W(-no:po) = Ws(-no:po)

      ! right boundary zone
      do k = po-no+1, po
        W(k) = W(k) + Ws(2*po + 1 - k)
      end do

    end if

  end associate

end subroutine GetSubdomainOperators

!-------------------------------------------------------------------------------
!> Finalization of a SchwarzOperator object

subroutine Delete_SchwarzOperator(this)
  type(SchwarzOperator), intent(inout) :: this !< Schwarz operator

  if (allocated(this%S1)) deallocate(this%S1)
  if (allocated(this%S2)) deallocate(this%S2)
  if (allocated(this%S3)) deallocate(this%S3)

  if (allocated(this%V1)) deallocate(this%V1)
  if (allocated(this%V2)) deallocate(this%V2)
  if (allocated(this%V3)) deallocate(this%V3)

  if (allocated(this%W1)) deallocate(this%W1)
  if (allocated(this%W2)) deallocate(this%W2)
  if (allocated(this%W3)) deallocate(this%W3)

end subroutine Delete_SchwarzOperator

!-------------------------------------------------------------------------------
!> Identification of element-boundary configurations

subroutine IdentifyConfiguration(element, bc, ec)
  type(MeshElement), intent(in) :: element(:)   !< mesh elements
  character, intent(in)  :: bc(:)               !< boundary conditions {'D','N'}
  integer,   intent(out) :: ec(3,size(element)) !< element configurations

  integer, parameter :: config(3,3) = reshape( [  1, &!  (1,1) = ' ', ' '
                                                  2, &!  (2,1) = 'D', ' '
                                                  4, &!  (3,1) = 'N', ' '
                                                  3, &!  (1,2) = ' ', 'D'
                                                  6, &!  (2,2) = 'D', 'D'
                                                  8, &!  (3,2) = 'N', 'D'
                                                  5, &!  (1,3) = ' ', 'N'
                                                  7, &!  (2,3) = 'D', 'N'
                                                  9  &!  (3,3) = 'N', 'N'
                                               ],    &
                                               [3,3] )

  character :: bc_face(size(bc))
  integer   :: e, f, c(6)

  do e = 1, size(element)

    where(element(e) % face % boundary > 0)
      bc_face = bc( element(e) % face % boundary )
    elsewhere
      bc_face = ''
    end where

    do f = 1, 6
      select case(bc_face(f))
      case('D')
        c(f) = 2
      case('N')
        c(f) = 3
      case default
        c(f) = 1
      end select
    end do

    ec(1,e) = config( c(1), c(2) )
    ec(2,e) = config( c(3), c(4) )
    ec(3,e) = config( c(5), c(6) )

  end do

end subroutine IdentifyConfiguration

!-------------------------------------------------------------------------------
!>

subroutine Iteration_C(this, mesh, lambda, nu, bc, u, f, i_max, r_red, r_max)

  ! arguments ..................................................................

  class(SchwarzOperator), intent(in) :: this  !< Schwarz operator
  class(MeshPartition),   intent(in) :: mesh  !< mesh partition

  real(RNP), intent(in)    :: lambda          !< Helmholtz parameter
  real(RNP), intent(in)    :: nu              !< constant diffusivity
  character, intent(in)    :: bc(:)           !< boundary conditions {P,D,N}
  real(RNP), intent(inout) :: u(:,:,:,:)      !< approximate solution
  real(RNP), intent(in)    :: f(:,:,:,:)      !< right hand side
  integer,   intent(in)    :: i_max           !< max number of iterations

  real(RNP), optional, intent(in) :: r_red    !< targeted residual reduction
  real(RNP), optional, intent(in) :: r_max    !< admissible residual

  ! local variables ............................................................

  !! save attribute serves to enforce sharing with OpenMP !!

  real(RNP), allocatable, save :: r (:,:,:,:)  ! residual, including ghosts
  real(RNP), allocatable, save :: u_s(:,:,:,:) ! solution of subsystems
  real(RNP), allocatable, save :: f_s(:,:,:,:) ! RHS of subsystems
  real(RNP), allocatable, save :: nu_s(:)      ! subdomain diffusivities
  integer,   allocatable, save :: ec(:,:)      ! element-boundary configurations

  type(ElementTransferBuffer), allocatable, save :: buf_r
  type(ElementTransferBuffer), allocatable, save :: buf_u_s

! procedure(TPO_Schwarz_Proc),     pointer :: SchwarzGenOP
  procedure(TPO_Schwarz_Iso_Proc), pointer :: SchwarzIsoOP

  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP) :: rr
  integer :: ne, ng, np(3), no(3), ns(3), nc
  integer :: i

  ! initialization .............................................................

  ! dimensions
  ne = mesh % ne              ! number of local elements
  ng = mesh % ng              ! number of ghost elements
  np = this % eop % po + 1    ! number of element points per direction
  no = this % no              ! overlapped point layers
  ns = np + 2*no              ! subdomain extensions
  nc = size(boundary_configuration, 2)

  ! workspace
  !$omp single
  allocate( r   (np(1), np(2), np(3), ne+ng) )
  allocate( u_s (ns(1), ns(2), ns(3), ne+ng) )
  allocate( f_s (ns(1), ns(2), ns(3), ne   ) )
  allocate( nu_s(ne), ec(3, ne) )
  allocate( buf_r )
  allocate( buf_u_s )
  !$omp end single

  call AssignScalar(r  , ZERO)
  call AssignScalar(u_s, ZERO)

  ! subdomain diffusivity -- constant, so far
  call AssignScalar(nu_s, nu)

  ! transfer buffers
  call buf_r   % New(mesh, r,  no)
  call buf_u_s % New(mesh, u_s, no)

  ! element-boundary configurations
  call IdentifyConfiguration(mesh%element, bc, ec)

  ! subdomain operator
  if (this % isotropic) then
    call TPO_Schwarz_Iso_Assign(ns(1), SchwarzIsoOP)
  else
!   call TPO_Schwarz_Assign(ns(1), ns(2), ns(3), SchwarzGenOP)
    call Error('Iteration_C', &
               'Anisotropic Schwarz operator not implemented yet', &
               'CART__DG_Elliptic_CI_Schwarz')
  end if

  ! termination condition
  !$omp single
  if (present(r_max)) then
    rr_term = max(ZERO, r_max)**2
  else
    rr_term = ZERO
  end if
  !$omp end single

  ! Schwarz iterations .........................................................

  do i = 1, i_max

    call EllipticResidual(mesh, this%eop, lambda, nu, bc, u, f, r(:,:,:,:ne))

    ! termination check
    if (present(r_red)) then
      rr = ScalarProduct(r(:,:,:,:ne), r(:,:,:,:ne), mesh%comm)
      !$omp master
      if (mesh%part == 0) then
        if (i == 1) then
          rr_term  = max(rr_term, max(ZERO, sqrt(rr) * r_red)**2)
        end if
        converged = rr <= rr_term
      end if
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier
    end if
    if (converged) exit

    call RestrictToSubdomains(mesh, this, buf_r, r, f_s)

    if (this % isotropic) then
      call SchwarzIsoOP( ns(1), nc, ne,                   &
                         this % S1,                       &
                         this % V1, this % V2, this % V3, &
                         this % W1,                       &
                         ec, mesh%dx, lambda,             &
                         nu_s, f_s, u_s                   )

    else
!     call SchwarzGenOP( ns(1), ns(2), ns(3), nc, ne,     &
!                        this % S1, this % S2, this % S3, &
!                        this % D1, this % D2, this % D3, &
!                        this % W1, this % W2, this % W3, &
!                        ec, f_s, u_s                       )
    end if

    call MergeFromSubdomains(mesh, this, buf_u_s, u_s, u)

  end do

  ! clean-up ...................................................................

  !$omp single
  deallocate(r, u_s, f_s, nu_s, ec)
  deallocate(buf_r, buf_u_s)
  !$omp end single

end subroutine Iteration_C

!===============================================================================
! Restriction and merging of subdomain data

!-------------------------------------------------------------------------------
!> Restrict a mesh variable to subdomain variable
!>
!> The strategy is as a follows:
!>
!>   1. Send `v` from linked local element to ghosts in remote partitions
!>   2. Assign `v` to subdomain core regions, i.e. 1:1 copy to `vs`
!>   3. Assign data from remote masters to ghosts
!>   4. Assign local and ghost data to overlap regions of `vs`
!>
!> The mesh variable must be dimensioned as `v(0:po,0:po,0:po,mesh%ne+mesh%ng)`.

subroutine RestrictToSubdomains(mesh, schwarz, buf_v, v, vs)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh    !< mesh partition
  class(SchwarzOperator),       intent(in)    :: schwarz !< Schwarz operator
  class(ElementTransferBuffer), intent(inout) :: buf_v   !< buffer for v

  real(RNP), intent(inout) :: v(0:,0:,0:,:)  !< extended mesh variable
  real(RNP), intent(out)   :: vs(:,:,:,:)    !< restricted variable

  ! local variables ............................................................

  integer :: po
  integer :: e, l, m, n, o
  integer :: i, ie0, ie1, is0, is1
  integer :: j, je0, je1, js0, js1
  integer :: k, ke0, ke1, ks0, ks1

  associate(no => schwarz % no)

    ! initialization ...........................................................

    po = ubound(v,1)

    call AssignScalar(vs, ZERO)

    ! offsets of subdomain point indices
    is0 = no(1) + 1;  is1 = po + is0
    js0 = no(2) + 1;  js1 = po + js0
    ks0 = no(3) + 1;  ks1 = po + ks0

    ! offsets of element point indices
    ie0 = -1;  ie1 = po + 1
    je0 = -1;  je1 = po + 1
    ke0 = -1;  ke1 = po + 1

    ! start transfer ...........................................................

    call buf_v % ToGhost_Transfer(mesh, v, tag=1000)

    ! assign core regions .................................................... ..

    do e = 1, mesh%ne
      do k = 0, po
      do j = 0, po
      do i = 0, po
        vs(is0 + i, js0 + j, ks0 + k, e) = v(i,j,k,e)
      end do
      end do
      end do
    end do

    ! assign remote data to ghosts .............................................

    ! for all ghosts: v = v[masters]
    call buf_v % ToGhost_Merge(v, alpha=ZERO, beta=ONE)

    ! fill overlap regions .....................................................

    ! offsets of subdomain point indices
    ! is0, js0, ks0 = 0 : not used
    ! is1, js1, ks1     : unchanged

    ! offsets of element point indices
    ie0 = -1;  ie1 = po - no(1)
    je0 = -1;  je1 = po - no(2)
    ke0 = -1;  ke1 = po - no(3)

    associate(element => mesh % element)

      do e = 1, mesh % ne

        ! face 1: -x1 <--> west
        o = element(e) % face(1) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            vs(l, js0 + j, ks0 + k, e) = v(ie1 + l, j, k, o)
          end do
          end do
          end do
        end if

        ! face 2:  x1 <--> east
        o = element(e) % face(2) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            vs(is1 + l, js0 + j, ks0 + k, e) = v(ie0 + l, j, k, o)
          end do
          end do
          end do
        end if

        ! face 3: -x2 <--> south
        o = element(e) % face(3) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, m, ks0 + k, e) = v(i, je1 + m, k, o)
          end do
          end do
          end do
        end if

        ! face 4:  x2 <--> north
        o = element(e) % face(4) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, js1 + m, ks0 + k, e) = v(i, je0 + m, k, o)
          end do
          end do
          end do
        end if

        ! face 5: -x3 <--> bottom
        o = element(e) % face(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            vs(is0 + i, js0 + j, n, e) = v(i, j, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! face 6:  x3 <--> top
        o = element(e) % face(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            vs(is0 + i, js0 + j, ks1 + n, e) = v(i, j, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! edge  1: (x1  , south, bottom)
        o = element(e) % edge( 1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, m, n, e) = v(i, je1 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  2: (x1  , north, bottom)
        o = element(e) % edge( 2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, js1 + m, n, e) = v(i, je0 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  3: (x1  , south, top   )
        o = element(e) % edge( 3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, m, ks1 + n, e) = v(i, je1 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! edge  4: (x1  , north, top   )
        o = element(e) % edge( 4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            vs(is0 + i, js1 + m, ks1 + n, e) = v(i, je0 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! edge  5: (west, x2   , bottom)
        o = element(e) % edge( 5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            vs(l, js0 + j, n, e) = v(ie1 + l, j, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  6: (east, x2   , bottom)
        o = element(e) % edge( 6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            vs(is1 + l, js0 + j, n, e) = v(ie0 + l, j, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  7: (west, x2   , top   )
        o = element(e) % edge( 7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            vs(l, js0 + j, ks1 + n, e) = v(ie1 + l, j, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! edge  8: (east, x2   , top   )
        o = element(e) % edge( 8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            vs(is1 + l, js0 + j, ks1 + n, e) = v(ie0 + l, j, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! edge  9: (west, south, x3    )
        o = element(e) % edge( 9) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, m, ks0 + k, e) = v(ie1 + l, je1 + m, k, o)
          end do
          end do
          end do
        end if

        ! edge 10: (east, south, x3    )
        o = element(e) % edge(10) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, m, ks0 + k, e) = v(ie0 + l, je1 + m, k, o)
          end do
          end do
          end do
        end if

        ! edge 11: (west, north, x3    )
        o = element(e) % edge(11) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, js1 + m, ks0 + k, e) = v(ie1 + l, je0 + m, k, o)
          end do
          end do
          end do
        end if

        ! edge 12: (east, north, x3    )
        o = element(e) % edge(12) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, js1 + m, ks0 + k, e) = v(ie0 + l, je0 + m, k, o)
          end do
          end do
          end do
        end if

        ! vertex 1: (west, south, bottom)
        o = element(e) % vertex(1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, m, n, e) = v(ie1 + l, je1 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 2: (east, south, bottom)
        o = element(e) % vertex(2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, m, n, e) = v(ie0 + l, je1 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 3: (west, north, bottom)
        o = element(e) % vertex(3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, js1 + m, n, e) = v(ie1 + l, je0 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 4: (east, north, bottom)
        o = element(e) % vertex(4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, js1 + m, n, e) = v(ie0 + l, je0 + m, ke1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 5: (west, south, top   )
        o = element(e) % vertex(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, m, ks1 + n, e) = v(ie1 + l, je1 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 6: (east, south, top   )
        o = element(e) % vertex(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, m, ks1 + n, e) = v(ie0 + l, je1 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 7: (west, north, top   )
        o = element(e) % vertex(7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(l, js1 + m, ks1 + n, e) = v(ie1 + l, je0 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 8: (east, north, top   )
        o = element(e) % vertex(8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            vs(is1 + l, js1 + m, ks1 + n, e) = v(ie0 + l, je0 + m, ke0 + n, o)
          end do
          end do
          end do
        end if

      end do
    end associate

    ! clean-up .................................................................

    call buf_v % ToGhost_Finish()

  end associate

end subroutine RestrictToSubdomains

!-------------------------------------------------------------------------------
!> Merge subdomain contributions into mesh variable
!>
!> The strategy is as a follows:
!>
!>   1. Apply weights to `vs`, taking into account the boundary conditions
!>   2. Send `vs` from masters to ghosts
!>   3. Add `vs` from subdomain core regions to `v`
!>   4. Assign data from remote masters to `vs` ghosts
!>   5. Merge `vs` from local and ghost overlap regions into `v`
!>
!> For transfer from masters to ghosts the subdomain variable must be shaped
!> `vs(ns(1), ns(2), ns(3), mesh%ne + mesh%ng)`.

subroutine MergeFromSubdomains(mesh, schwarz, buf_vs, vs, v)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh    !< mesh partition
  class(SchwarzOperator),       intent(in)    :: schwarz !< Schwarz operator
  class(ElementTransferBuffer), intent(inout) :: buf_vs  !< buffer for vs

  real(RNP), intent(inout) :: vs(:,:,:,:)    !< subdomain variables
  real(RNP), intent(inout) :: v(0:,0:,0:,:)  !< mesh variable

  ! local variables ............................................................

  integer :: po, no(3), ns(3)
  integer :: e, l, m, n, o
  integer :: i, ie0, ie1, is0, is1
  integer :: j, je0, je1, js0, js1
  integer :: k, ke0, ke1, ks0, ks1

  ! initialization ...........................................................

  ! dimensions
  po = ubound(v,1)
  no = schwarz % no
  ns = po + 1 + 2 * no

  ! offsets of subdomain point indices
  is0 = no(1) + 1;  is1 = ns(1) - no(1)
  js0 = no(2) + 1;  js1 = ns(2) - no(2)
  ks0 = no(3) + 1;  ks1 = ns(3) - no(3)

  ! offsets of element point indices
  ie0 = -1;  ie1 = po - no(1)
  je0 = -1;  je1 = po - no(2)
  ke0 = -1;  ke1 = po - no(3)

  ! send vs to ghosts ........................................................

  call buf_vs % ToGhost_Transfer(mesh, vs, tag=2000)

  ! add vs core regions to v .................................................

  do e = 1, mesh%ne
    do k = 0, po
    do j = 0, po
    do i = 0, po
      v(i,j,k,e) = v(i,j,k,e) + vs(is0 + i, js0 + j, ks0 + k, e)
    end do
    end do
    end do
  end do

  ! copy received data to ghosts .............................................

  call buf_vs % ToGhost_Merge(vs, alpha=ZERO, beta=ONE)

  ! merge vs from overlap regions into v .....................................

    associate(element => mesh % element)

      do e = 1, mesh % ne

        ! face 1: -x1 <--> west
        o = element(e) % face(1) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, k, e) = v (ie0 + l,       j,       k, e) &
                                + vs(is1 + l, js0 + j, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 2:  x1 <--> east
        o = element(e) % face(2) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, k, e) = v (ie1 + l,       j,       k, e) &
                                + vs(      l, js0 + j, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 3: -x2 <--> south
        o = element(e) % face(3) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, k, e) = v (      i, je0 + m,       k, e) &
                                + vs(is0 + i, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 4:  x2 <--> north
        o = element(e) % face(4) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, k, e) = v (      i, je1 + m,       k, e) &
                                + vs(is0 + i,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 5: -x3 <--> bottom
        o = element(e) % face(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            v(i, j, ke0 + n, e) = v (      i,       j, ke0 + n, e) &
                                + vs(is0 + i, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! face 6:  x3 <--> top
        o = element(e) % face(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            v(i, j, ke1 + n, e) = v (        i,         j, ke1 + n, e) &
                                + vs(is0 + i, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  1: (x1  , south, bottom)
        o = element(e) % edge( 1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, ke0 + n, e) = v (      i, je0 + m, ke0 + n, e) &
                                      + vs(is0 + i, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  2: (x1  , north, bottom)
        o = element(e) % edge( 2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, ke0 + n, e) = v (      i, je1 + m, ke0 + n, e) &
                                      + vs(is0 + i,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  3: (x1  , south, top   )
        o = element(e) % edge( 3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, ke1 + n, e) = v (      i, je0 + m, ke1 + n, e) &
                                      + vs(is0 + i, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! edge  4: (x1  , north, top   )
        o = element(e) % edge( 4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, ke1 + n, e) = v (      i, je1 + m, ke1 + n, e) &
                                      + vs(is0 + i,       m,       n, o)
          end do
          end do
          end do
        end if

        ! edge  5: (west, x2   , bottom)
        o = element(e) % edge( 5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, ke0 + n, e) = v (ie0 + l,       j, ke0 + n, e) &
                                      + vs(is1 + l, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  6: (east, x2   , bottom)
        o = element(e) % edge( 6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, ke0 + n, e) = v (ie1 + l,       j, ke0 + n, e) &
                                      + vs(      l, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  7: (west, x2   , top   )
        o = element(e) % edge( 7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, ke1 + n, e) = v (ie0 + l,       j, ke1 + n, e) &
                                      + vs(is1 + l, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  8: (east, x2   , top   )
        o = element(e) % edge( 8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, ke1 + n, e) = v (ie1 + l,       j, ke1 + n, e) &
                                      + vs(      l, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  9: (west, south, x3    )
        o = element(e) % edge( 9) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, k, e) = v (ie0 + l, je0 + m,       k, e) &
                                      + vs(is1 + l, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 10: (east, south, x3    )
        o = element(e) % edge(10) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, k, e) = v (ie1 + l, je0 + m,       k, e) &
                                      + vs(      l, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 11: (west, north, x3    )
        o = element(e) % edge(11) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, k, e) = v (ie0 + l, je1 + m,       k, e) &
                                      + vs(is1 + l,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 12: (east, north, x3    )
        o = element(e) % edge(12) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, k, e) = v (ie1 + l, je1 + m,       k, e) &
                                      + vs(      l,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! vertex 1: (west, south, bottom)
        o = element(e) % vertex(1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, ke0 + n, e) = v (ie0 + l, je0 + m, ke0 + n, e) &
                                            + vs(is1 + l, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 2: (east, south, bottom)
        o = element(e) % vertex(2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, ke0 + n, e) = v (ie1 + l, je0 + m, ke0 + n, e) &
                                            + vs(      l, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 3: (west, north, bottom)
        o = element(e) % vertex(3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, ke0 + n, e) = v (ie0 + l, je1 + m, ke0 + n, e) &
                                            + vs(is1 + l,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 4: (east, north, bottom)
        o = element(e) % vertex(4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, ke0 + n, e) = v (ie1 + l, je1 + m, ke0 + n, e) &
                                            + vs(      l,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 5: (west, south, top   )
        o = element(e) % vertex(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, ke1 + n, e) = v (ie0 + l, je0 + m, ke1 + n, e) &
                                            + vs(is1 + l, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 6: (east, south, top   )
        o = element(e) % vertex(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, ke1 + n, e) = v (ie1 + l, je0 + m, ke1 + n, e) &
                                            + vs(      l, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 7: (west, north, top   )
        o = element(e) % vertex(7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, ke1 + n, e) = v (ie0 + l, je1 + m, ke1 + n, e) &
                                            + vs(is1 + l,       m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 8: (east, north, top   )
        o = element(e) % vertex(8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, ke1 + n, e) = v (ie1 + l, je1 + m, ke1 + n, e) &
                                            + vs(      l,       m,       n, o)
          end do
          end do
          end do
        end if

      end do
    end associate

  ! clean-up .................................................................

  call buf_vs % ToGhost_Finish()

end subroutine MergeFromSubdomains

!===============================================================================

end module CART__DG_Elliptic_CI_Schwarz
