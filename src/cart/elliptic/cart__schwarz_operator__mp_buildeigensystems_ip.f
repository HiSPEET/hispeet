!> summary:  Build 1D Schwarz eigensystems for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2018/11/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build 1D Schwarz eigensystems for IP/DG-SEM
!===============================================================================

submodule(CART__Schwarz_Operator) MP_BuildEigensystems_IP
  use Constants,     only: ZERO, ONE, HALF
  use Eigenproblems, only: SolveGeneralizedEigenproblem
  use Schwarz_Weighting
  use IP_Element_Operators_1D
  implicit none

contains

!-------------------------------------------------------------------------------
!> Build 1D Schwarz eigensystems for IP/DG-SEM

module subroutine BuildEigensystems_IP(this, eop, delta, no_min, weighting)

  ! arguments ..................................................................

  class(SchwarzOperator3D),     intent(inout) :: this !< Schwarz operator
  class(IP_ElementOperators1D), intent(in)    :: eop  !< IP-DG SE operators

  real(RNP),         intent(in) :: delta(3)  !< relative overlap
  integer, optional, intent(in) :: no_min    !< min overlap in points [-1]
  integer, optional, intent(in) :: weighting !< weighting method      [ 5]

  ! local variables ............................................................

  real(RNP), allocatable :: W1(:), W2(:), W3(:)

  character :: bc(2)
  integer :: i, j, k

  associate(po => eop%po, Ms => eop%w)

    ! parameters ...............................................................

    do i = 1, 3
      this % no(i) = count(eop%x <= 2*delta(i) - 1)
    end do

    if (present(no_min)) then
      this % no = max(this%no, min(no_min, po+1))
    end if

    this % isotropic = all(this%no(2:3) == this%no(1))

    this % n1 = po + 1 + 2 * this%no(1)
    this % n2 = po + 1 + 2 * this%no(2)
    this % n3 = po + 1 + 2 * this%no(3)

    this % nc = num_boundary_types ** 2

    ! preparations .............................................................

    associate(n1 => this%n1, n2 => this%n2, n3 => this%n3, nc => this%nc)

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

    end associate

    ! eigenvectors, eigenvalues and weights ....................................

    do j = 1, num_boundary_types

      bc(2) = boundary_type(j)

      do i = 1, num_boundary_types

        bc(1) = boundary_type(i)
        k = i + num_boundary_types * (j - 1)

        call GetSubdomainOperators( eop,                 &
                                    no = this%no(1),     &
                                    bc = bc,             &
                                    Ws = W1,             &
                                    S  = this%S1(:,:,k), &
                                    V  = this%V1(:,k),   &
                                    W  = this%W1(:,k)    )

        if (this%isotropic) then

          this % S2 = this % S1
          this % V2 = this % V1
          this % W2 = this % W1

          this % S3 = this % S1
          this % V3 = this % V1
          this % W3 = this % W1

        else

          call GetSubdomainOperators( eop,                 &
                                      no = this%no(2),     &
                                      bc = bc,             &
                                      Ws = W2,             &
                                      S  = this%S2(:,:,k), &
                                      V  = this%V2(:,k),   &
                                      W  = this%W2(:,k)    )

          call GetSubdomainOperators( eop,                 &
                                      no = this%no(3),     &
                                      bc = bc,             &
                                      Ws = W3,             &
                                      S  = this%S3(:,:,k), &
                                      V  = this%V3(:,k),   &
                                      W  = this%W3(:,k)    )
        end if

      end do
    end do

  end associate

end subroutine BuildEigensystems_IP

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
!> Although no exterior node layers exist at boundaries, the corresponding
!> entries are retained for regularity. They are, however, set to `0` in the
!> eigenvectors and weights. The corresponding eigenvalues they are set `1` in
!> order to avoid floating points exceptions when used as a divisor.

subroutine GetSubdomainOperators(eop, no, bc, Ws, S, V, W)
  type(IP_ElementOperators1D), intent(in) :: eop !< IP-DG element operators
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

  ! stiffness matrix for interior element, assuming dx=1
  call eop % GetStiffnessMatrix(dx, bc=ii, Le=Le_ii)

  ! stiffness matrix for given boundary conditions, assuming dx=1
  call eop % GetStiffnessMatrix(dx, bc=bc, Le=Le_bc)

  ! initialization: eigenvalues V set to 1 to avoid division by zero in
  ! in case of truncated overlap zones
  S = 0
  V = 1
  W = 0

  associate(Ms => eop % W)

    if (all(bc == ' ')) then

      ! interior-interior-interior configuration -------------------------------

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

      ! boundary-interior-boundary configuration -------------------------------

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

      ! boundary-interior-interior configuration -------------------------------

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

      ! interior-interior-boundary configuration -------------------------------

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

!===============================================================================

end submodule MP_BuildEigensystems_IP
