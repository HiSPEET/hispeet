!> summary:  Schwarz operator for elliptic equations
!> author:   Joerg Stiller
!> date:     2022/02/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__3D
  use Kind_Parameters
  use Constants, only: ONE, ZERO

  !-----------------------------------------------------------------------------
  !> Schwarz operator
  !>
  !> <add details>

  type DG_SchwarzOperator_3D

    integer :: po = -1     !< polynomial order
    integer :: no = -1     !< number of overlapped layers
    integer :: nc = -1     !< number of 1D configurations
    logical :: restrictive !< T/F for ex/including neighbor results

    real(RNP), allocatable :: S(:,:,:) !< eigenvectors per config
    real(RNP), allocatable :: V(:,:)   !< eigenvalues  per config
    real(RNP), allocatable :: W(:,:)   !< weights      per config

  end type DG_SchwarzOperator_3D

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the Schwarz operator
  !>
  !> <add details>

  subroutine Init_DG_SchwarzOperator_3D(this, weighting, xi, bc, L, M)
    class(DG_SchwarzOperator_3D), intent(inout) :: this
    integer,   intent(in) :: weighting !< weighting method {0,1,3,5,7,9}
    real(RNP), intent(in) :: xi(0:)    !< standard element points
    logical,   intent(in) :: bc(:,:)   !< boundary configurations
    real(RNP), intent(in) :: L(:,:,:)  !< standard subdomain operator matrices
    real(RNP), intent(in) :: M(:,:)    !< standard subdomain mass matrices

    ! local variables ..........................................................

    real(RDP), allocatable :: L0(:,:), M0(:), S0(:,:), V0(:), W0(:) ! interior
    real(RDP), allocatable :: L1(:,:), M1(:), S1(:,:), V1(:)        ! one boundary
    real(RDP), allocatable :: L2(:,:), M2(:), S2(:,:), V2(:)       ! two boundaries

    integer :: nc, no, np, ns, po
    integer :: n0, n1, n2
    integer :: c, j, k

    ! dimensions ...............................................................

    po = ubound(xi,1)   ! polynomial order
    np = po + 1         ! number of element points per direction
    ns = size(M,1)      ! number of subdomain points per direction
    nc = size(M,2)      ! number of configurations
    no = (ns - np) / 2  ! number of overlapped points

    n0 = ns             ! number of points for interior-interior configuration
    n1 = np + no        ! number of points for interior-boundary configuration
    n2 = np             ! number of points for boundary-boundary configuration

    ! parameters ...............................................................

    this % po = po
    this % no = no
    this % nc = nc
    this % restrictive = weighting == 9 ! top-hat weighting !

    ! eigensystems .............................................................

    if (allocated(this % S)) deallocate(this % S)
    if (allocated(this % V)) deallocate(this % V)
    if (allocated(this % W)) deallocate(this % W)

    allocate(this % S(ns,ns,nc), source = ZERO)
    allocate(this % V(ns,nc)   , source = ONE )
    allocate(this % W(ns,nc)   , source = ZERO)

    allocate( L0(n0,n0), M0(n0), S0(n0,n0), V0(n0), W0(n0) )
    allocate( L1(n1,n1), M1(n1), S1(n1,n1), V1(n1) )
    allocate( L2(n2,n2), M2(n2), S2(n2,n2), V2(n2) )

    call WeightDistribution(xi, no, weighting, W0)

    do c = 1, nc

      if (bc(1,c) and bc(2,c)) then ! b-b  ↔︎  subdomain truncated on both sides

        j = no + 1
        k = no + np

        L2 = L(j:k, j:k, c)
        M2 = M(j:k, c)

        call SolveGeneralizedEigenproblem(L2, M2, V2, S2)

        this % S(j:k, j:k, c) = S2
        this % V(j:k, c)      = V2
        this % W(j:k, c)      = W0

        ! adjust weights at left boundary
        do i = 0, no-1
          this % W(j+i,c) = this % W(j+i,c) + W0(j-i-1)
        end do

        ! adjust weights at right boundary
        do i = 0, no-1
          this % W(k-i,c) = this % W(k-i,c) + W0(k+i+1)
        end do

      else if (bc(1,c)) then ! b-i  ↔︎  subdomain truncated on the left side

        j = no + 1
        k = no + np + no

        L1 = L(j:k, j:k, c)
        M1 = M(j:k, c)

        call SolveGeneralizedEigenproblem(L1, M1, V1, S1)

        this % S(j:k, j:k, c) = S1
        this % V(j:k, c)      = V1
        this % W(j:k, c)      = W0

        ! adjust weights at left boundary
        do i = 0, no-1
          this % W(j+i,c) = this % W(j+i,c) + W0(j-i-1)
        end do

      else if (bc(2,c)) then ! i-b  ↔︎  subdomain truncated on the right side

        j =      1
        k = no + np

        L1 = L(j:k, j:k, c)
        M1 = M(j:k, c)

        call SolveGeneralizedEigenproblem(L1, M1, V1, S1)

        this % S(j:k, j:k, c) = S1
        this % V(j:k, c)      = V1
        this % W(j:k, c)      = W0

        ! adjust weights at right boundary
        do i = 0, no-1
          this % W(k-i,c) = this % W(k-i,c) + W0(k+i+1)
        end do

      else ! i-i  ↔  interior subdomain wi︎th overlap on both sides

        j =      1
        k = no + np + no

        L0 = L(j:k, j:k, c)
        M0 = M(j:k, c)

        call SolveGeneralizedEigenproblem(L0, M0, V0, S0)

        this % S(j:k, j:k, c) = S0
        this % V(j:k, c)      = V0
        this % W(j:k, c)      = W0

      end if

    end do

  end subroutine Init_DG_SchwarzOperator_3D

  !==============================================================================

end module DG__Schwarz_Operator__3D
