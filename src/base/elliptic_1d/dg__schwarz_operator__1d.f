!> summary:  DG Schwarz operator for elliptic equations
!> author:   Joerg Stiller
!> date:     2023/03/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__1D
  use Kind_Parameters
  use Constants
  use Schwarz_Weighting
  use DG__Element_Operators__1D

  implicit none
  private

  public :: DG_SCHWARZ_BC_1D
  public :: DG_SchwarzOptions_1D
  public :: DG_SchwarzOperator_1D

  !-----------------------------------------------------------------------------
  !> Supported boundary conditions
  !>
  !> Periodic boundaries are treated as interior: 'P' → ' '

  character, parameter :: DG_SCHWARZ_BC_1D(3) = [ ' ', 'D', 'N']

  !-----------------------------------------------------------------------------
  !> Options for initializing the Schwarz operator

  type DG_SchwarzOptions_1D
    integer :: no = 1         !< number of overlapped points
    integer :: weighting = 5  !< weighting method {0,1,3,5,7,9}
  end type DG_SchwarzOptions_1D

  !-----------------------------------------------------------------------------
  !> Schwarz operator
  !>
  !> In the Schwarz method we consider a subdomain comprising one element
  !> located in its center and `no` layers of collocation points which are
  !> are adopted from the adjoining elements.
  !>
  !> The Schwarz operator is the inverse of the truncated diffusion operator,
  !> which is given in tensor-product form by
  !>
  !>     A  =  λ ∆x M  +  ν ∆x⁻¹ L
  !>
  !> where `M` is the 1D mass matrix, `L` the corresponding stiffness matrix,
  !> λ the Helmholtz parameter and ν the diffusivity coefficient.
  !> The operators `M` and `L` are normalized to unit mesh spacing and depend
  !> on the following parameters
  !>
  !>   -  polynomial order
  !>   –  further discretization parameters such as penalty factors
  !>   -  number of overlapped points `no`
  !>   -  boundary conditions.
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
  !> solved, yielding the matrix of right eigenvectors `S` and the diagonal
  !> matrix of eigenvalues `Λ` such that
  !>
  !>     Sᵀ L S = Λ
  !>     Sᵀ M S = I
  !>
  !> where `I` ist the matching unit matrix. The inverse Helmholtz operator
  !> can be then expressed in the tensor-product form
  !>
  !>     A⁻¹  =  S D⁻¹ Sᵀ
  !>
  !> with the diagonal matrix
  !>
  !>     D  =  λ ∆x I  +  ν ∆x⁻¹ Λ
  !>
  !>
  !> As λ and ν can vary, `D⁻¹` is not stored, but computed on the fly using
  !> the 1D eigenvalues, which are kept within the operator.
  !>
  !> Before assembling the global correction to an approximate solution, the
  !> subdomain correction is weighted according to
  !>
  !>     Δu = W (A⁻¹ r)
  !>
  !> where `W` is the diagonal matrix containing the weights, depending on
  !> the element-boundary configuration.

  type DG_SchwarzOperator_1D

    integer :: po = -1     !< polynomial order
    integer :: no = -1     !< number of overlapped layers
    integer :: nc = -1     !< number of 1D configurations
    logical :: restrictive !< T/F for ex/including neighbor results

    real(RNP), allocatable :: S(:,:,:) !< eigenvectors per config
    real(RNP), allocatable :: V(:,:)   !< eigenvalues  per config
    real(RNP), allocatable :: W(:,:)   !< weights      per config

  contains

    procedure :: ConfigureSubdomains
    procedure :: RestrictResidual
    procedure :: MergeCorrections
    generic   :: Apply => Apply_R
    procedure, private :: Apply_R

  end type DG_SchwarzOperator_1D

  ! constructor interface
  interface DG_SchwarzOperator_1D
    module procedure New_DG_SchwarzOperator_1D
  end interface

contains

  !=============================================================================
  ! Type-bound procedures of DG_SchwarzOperator_1D

  !-----------------------------------------------------------------------------
  !> Constructor

  function New_DG_SchwarzOperator_1D(opt, eop, r_nu_s) result(this)

    class(DG_SchwarzOptions_1D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    real(RNP),           optional, intent(in) :: r_nu_s !< ratio νˢ/(ν + νˢ) [0]

    type(DG_SchwarzOperator_1D) :: this

    call InitSchwarzOperator(this, opt, eop, r_nu_s)

  end function New_DG_SchwarzOperator_1D

  !-----------------------------------------------------------------------------
  !> Build the 1D eigenvalues, eigenvectors and weights

  subroutine InitSchwarzOperator(this, opt, eop, r_nu_s)
    class(DG_SchwarzOperator_1D),  intent(inout) :: this
    class(DG_SchwarzOptions_1D),   intent(in) :: opt    !< operator options
    class(DG_ElementOperators_1D), intent(in) :: eop    !< DG element operators
    real(RNP),           optional, intent(in) :: r_nu_s !< ratio νˢ/(ν + νˢ) [0]

    ! local variables ..........................................................

    real(RNP), allocatable :: Ws(:)
    character :: bc(2)

    integer :: nb, nc, no, ns, po
    integer :: i, j, k

    ! initialization ...........................................................

    po = eop % po                  ! polynomial order of elements
    nb = size(DG_SCHWARZ_BC_1D)    ! number of supported boundary conditions
    nc = nb ** 2                   ! number of 1D boundary configurations
    no = max(0, min(opt%no, po+1)) ! number of overlapped points
    ns = po + 1 + 2*no             ! number of subdomain points per direction

    allocate(Ws(ns))

    !$omp master

    this % po = po
    this % no = no
    this % nc = nc
    this % restrictive = opt % weighting == 9 .or. no == 0

    if (allocated(this % S  )) deallocate(this % S  )
    if (allocated(this % V  )) deallocate(this % V  )
    if (allocated(this % W  )) deallocate(this % W  )

    allocate(this % S(ns,ns,nc))
    allocate(this % V(ns,nc)   )
    allocate(this % W(ns,nc)   )

    !$omp end master
    !$omp barrier

    ! eigensystems and weights .................................................

    call WeightDistribution(eop%x, this%no, opt%weighting, Ws)

    !$omp do collapse(2)
    do j = 1, nb
    do i = 1, nb

      k = i + nb * (j - 1)

      bc(1) = DG_SCHWARZ_BC_1D(i)
      bc(2) = DG_SCHWARZ_BC_1D(j)

      call eop % Get_SchwarzSuboperators( this%no, bc, Ws, this%S(:,:,k)   &
                                        , this%V(:,k), this%W(:,k), r_nu_s )
    end do
    end do

  end subroutine InitSchwarzOperator

  !-----------------------------------------------------------------------------
  !> Configure subdomains for automatic operator selection

  subroutine ConfigureSubdomains(this, bc, cfg, mask)
    class(DG_SchwarzOperator_1D), intent(in) :: this
    character,         intent(in)  :: bc(2)  !< BC {'D','N','P'}
    integer,           intent(out) :: cfg(:) !< subdomain configurations
    logical, optional, intent(in)  :: mask(:)

    character :: bc_e(2)
    integer   :: e, ne

    ne = size(cfg)

    if (ne == 1) then
      cfg(1) = ConfigurationID(bc)

    else if (present(mask)) then
      ! masked configuration
      do e = 1, ne

        if (mask(e)) then

          ! subdomain bc on left side
          if (e == 1) then
            ! left boundary
            bc_e(1) = bc(1)
          else if (mask(e-1)) then
            ! left neighbor is active
            bc_e(1) = ' '
          else
            ! left neighbor is frozen
            bc_e(1) = 'D'
          end if

          ! subdomain bc on right side
          if (e == ne) then
            ! right boundary
            bc_e(2) = bc(2)
          else if (mask(e+1)) then
            ! right neighbor is active
            bc_e(2) = ' '
          else
            ! right neighbor is frozen
            bc_e(2) = 'D'
          end if

          cfg(e) = ConfigurationID(bc_e)

        else

          ! skip configuration of subdomain for frozen elements
          cfg(e) = 0

        end if
      end do

    else
      ! default configuration
      bc_e = [ bc(1),   ' ' ];  cfg(1     ) = ConfigurationID(bc_e)
      bc_e = [   ' ',   ' ' ];  cfg(2:ne-1) = ConfigurationID(bc_e)
      bc_e = [   ' ', bc(2) ];  cfg(  ne  ) = ConfigurationID(bc_e)

    end if

    ! silence compiler warnings ;)
    if (this % po > 0) return

  end subroutine ConfigureSubdomains

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

    cfg = i1 + (i2 - 1) * size(DG_SCHWARZ_BC_1D)

  end function ConfigurationID

  !-----------------------------------------------------------------------------
  !> Restrict mesh variable to subdomains

  subroutine RestrictResidual(this, periodic, r, rs)
    class(DG_SchwarzOperator_1D), intent(in) :: this
    logical,               intent(in)  :: periodic !< T/F for periodicity
    real(RNP), contiguous, intent(in)  :: r(0:,:)  !< mesh variable
    real(RDP), contiguous, intent(out) :: rs(:,:)  !< subdomain variable

    integer :: e, ne, np, ns

    associate(po => this%po, no => this%no)

      np = size(r ,1)
      ne = size(r ,2)
      ns = size(rs,1)

      ! core regions ...........................................................

      !$omp do
      do e = 1, ne
        rs(no+1:no+np, e) = r(:, e)
      end do
      !$omp end do nowait

      ! overlap regions ........................................................

      !$omp do
      do e = 2, ne-1
        rs(      1:no, e) = r(np-no:po  , e-1)
        rs(ns+1-no:ns, e) = r(    0:no-1, e+1)
      end do

      !$omp master
      if (no > 0) then
        if (ne > 1) then
          rs(ns+1-no:ns,  1) = r(    0:no-1,    2)
          rs(      1:no, ne) = r(np-no:po  , ne-1)
        end if
        if (periodic) then
          rs(      1:no,  1) = r(np-no:po  , ne)
          rs(ns+1-no:ns, ne) = r(    0:no-1,  1)
        else
          rs(      1:no,  1) = 0
          rs(ns+1-no:ns, ne) = 0
        end if
      end if
      !$omp end master

    end associate

  end subroutine RestrictResidual

  !-----------------------------------------------------------------------------
  !> Merge subdomain contributions into mesh variable

  subroutine MergeCorrections(this, periodic, us, u)
    class(DG_SchwarzOperator_1D), intent(in) :: this
    logical,               intent(in)    :: periodic !< T/F for periodicity
    real(RNP), contiguous, intent(inout) :: u(0:,:)  !< subdomain solutions
    real(RNP), contiguous, intent(in)    :: us(:,:)  !< mesh variable

    integer :: e, ne, np, ns

    associate(po => this%po, no => this%no)

      np = size(u ,1)
      ne = size(u ,2)
      ns = size(us,1)

      ! core regions ...........................................................

      !$omp do
      do e = 1, ne
        u(:, e) = u(:, e) + us(no+1:no+np, e)
      end do

      ! overlap regions ........................................................

      !$omp do
      do e = 2, ne-1
        u(    0:no-1, e) = u(    0:no-1, e) + us(ns+1-no:ns, e-1)
        u(np-no:po  , e) = u(np-no:po  , e) + us(      1:no, e+1)
      end do

      !$omp master
      if (no > 0) then
        if (ne > 1) then
          u(np-no:po  ,  1) = u(np-no:po  ,  1) + us(      1:no,    2)
          u(    0:no-1, ne) = u(    0:no-1, ne) + us(ns+1-no:ns, ne-1)
        end if
        if (periodic) then
          u(    0:no-1,  1) = u(    0:no-1,  1) + us(ns+1-no:ns, ne)
          u(np-no:po  , ne) = u(np-no:po  , ne) + us(      1:no,  1)
        end if
      end if
      !$omp end master

    end associate

  end subroutine MergeCorrections

  !-----------------------------------------------------------------------------
  !> Application with regular mesh
  !>
  !> For flexibility, ν is allowed to vary from subdomain to subdomain.

  subroutine Apply_R(this, cfg, dx, lambda, nu, rs, us)
    class(DG_SchwarzOperator_1D), intent(in) :: this
    integer,               intent(in)  :: cfg(:)  !< subdomain configurations
    real(RNP),             intent(in)  :: dx      !< ∆xᵉ
    real(RNP),             intent(in)  :: lambda  !< λ
    real(RNP),             intent(in)  :: nu(:)   !< ν = νᵖ + νˢ
    real(RNP), contiguous, intent(in)  :: rs(:,:) !< operand
    real(RNP), contiguous, intent(out) :: us(:,:) !< result

    real(RNP), parameter :: eps = epsilon(lambda)
    real(RNP) :: a0, a1, d, tmp, z(size(rs,1))
    integer   :: c, i, j, l, ne, ns

    associate(S => this%S, V => this%V, W => this%W)

      ns = size(rs,1)
      ne = size(rs,2)
      a0 = lambda * dx

      !$omp do
      do l = 1, ne

        if (cfg(l) == 0) then
          ! skip subdomains surrounding frozen elements
          us(:,l) = 0
        else

          c = cfg(l)
          a1 = nu(l) / dx

          do i = 1, ns
            tmp = 0
            do j = 1, ns
              tmp = tmp + S(j,i,c) * rs(j,l)
            end do
            d = a0 + a1 * V(i,c)
            z(i) = tmp / sign(max(abs(d),eps), d)
          end do

          do i = 1, ns
            tmp = 0
            do j = 1, ns
              tmp = tmp + S(i,j,c) * z(j)
            end do
            us(i,l) = W(i,c) * tmp
          end do

        end if
      end do

    end associate

  end subroutine Apply_R

  !=============================================================================

end module DG__Schwarz_Operator__1D
