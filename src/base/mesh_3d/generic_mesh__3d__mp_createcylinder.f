!> summary:  Creation of a 3d generic mesh in a cylinder
!> author:   Joerg Stiller
!> date:     2021/01/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!> @todo
!>   - flexible choice of radius and length
!>   – nonuniform element spacing
!==============================================================================

submodule(Generic_Mesh__3D) MP_CreateCylinder
  use Constants
  use Gauss_Jacobi
  implicit none

  ! curvature radius of one side of the inner quadrangle
  real(RHP), parameter :: RI = 1.1380711874576983496005629080698994_RHP

  ! distance of the midpoint of the osculating circle from origin
  real(RHP), parameter :: DI = RI - THIRD

contains

  !-----------------------------------------------------------------------------
  !> Create a uniform generic mesh for a cylinder of radius r and height h
  !>
  !> Restriction: nz = 1 or nz ≥ 3 in the periodic case.

  module subroutine CreateCylinder(mesh, r, h, nr, nz, po, periodic)
    class(GenericMesh_3D), intent(out) :: mesh  !< cylindrical 3d mesh
    real(RNP), intent(in) :: r   !< cylinder radius
    real(RNP), intent(in) :: h   !< cylinder height (axial extension)
    integer,   intent(in) :: nr  !< num intervals in radial section
    integer,   intent(in) :: nz  !< num intervals in axial  direction
    integer,   intent(in) :: po  !< polynomial order of mesh elements
    logical,   intent(in) :: periodic !< switch for axial periodicity

    integer   :: nv, ne, nb, np
    integer   :: i, j, k, l, d, f, m, n, o
    real(RNP) :: xc(0:po), cr, cz
    real(RNP) :: xi, eta, zeta

    ! prerequisites ............................................................

    ! number of vertices and elements and boundaries
    nv = ((nr+1)**2 + 4*nr**2) * (nz+1)
    ne = 5*nr**2 * nz
    nb = 3

    ! collocation points
    xc = LobattoPoints(po)
    np = size(xc)

    ! scaling factors
    cr = ONE / nr
    cz = ONE / nz

    ! allocate mesh components
    allocate(mesh%vertex(nv), mesh%element(ne), mesh%boundary(nb))

    ! set lexical numbering
    mesh % numbering = LEXICAL_NUMBERING

    ! vertices .................................................................

    do k = 0, nz
      zeta = TWO*k/nz - 1
      ! subdomain 1: inner quadrangle
      do concurrent (i=0:nr, j=0:nr)
        mesh%vertex(LinearVertexIndex(1, i, j, k, nr))%x = &
            PointCoordinates(1, 2*cr*i-1, 2*cr*j-1, zeta, r, h)
      end do
      ! subdomains 2-5: ring sections
      do concurrent (i=1:nr, j=1:nr, d=2:5)
        mesh%vertex(LinearVertexIndex(d, i, j, k, nr))%x = &
            PointCoordinates(d, 2*cr*i-1, 2*cr*j-1, zeta, r, h)
      end do
    end do

    ! elements .................................................................

    mesh%element%typ   = HEXAHEDRAL_ELEMENT
    mesh%element%basis = GAUSS_LOBATTO_BASIS
    mesh%element%order = po

    do k = 1, nz
    do d = 1, 5
    do j = 1, nr
    do i = 1, nr
      l = LinearElementIndex(d, i, j, k, nr)

      mesh%element(l)%id = l

      ! vertices
      mesh%element(l)%vertex = [ LinearVertexIndex(d, i-1, j-1, k-1, nr), &
                                 LinearVertexIndex(d, i  , j-1, k-1, nr), &
                                 LinearVertexIndex(d, i-1, j  , k-1, nr), &
                                 LinearVertexIndex(d, i  , j  , k-1, nr), &
                                 LinearVertexIndex(d, i-1, j-1, k  , nr), &
                                 LinearVertexIndex(d, i  , j-1, k  , nr), &
                                 LinearVertexIndex(d, i-1, j  , k  , nr), &
                                 LinearVertexIndex(d, i  , j  , k  , nr) ]

      ! points
      allocate(mesh%element(l)%x((po+1)**3, 3))
      xi   = cr * (2*i - 1)  -  1
      eta  = cr * (2*j - 1)  -  1
      zeta = cz * (2*k - 1)  -  1
      do o = 0, po
      do n = 0, po
      do m = 0, po
        mesh%element(l)%x(1 + m + np*(n + np*o),:) = &
            PointCoordinates(d, xi+cr*xc(m), eta+cr*xc(n), zeta+cz*xc(o), r, h)
      end do
      end do
      end do
    end do
    end do
    end do
    end do

    ! boundaries ...............................................................

    ! shell (lateral surface)
    mesh%boundary(1)%id = 1
    mesh%boundary(1)%name = 'shell'
    allocate(mesh%boundary(1)%face(4*nr*nz))
    f = 0
    i = nr
    do k = 1, nz
    do d = 2, 5
    do j = 1, nr
      f = f + 1
      mesh%boundary(1)%face(f)%element_id   = LinearElementIndex(d, i, j, k, nr)
      mesh%boundary(1)%face(f)%element_face = 2
    end do
    end do
    end do

    ! bottom
    mesh%boundary(2)%id = 2
    mesh%boundary(2)%name = 'bottom'
    allocate(mesh%boundary(2)%face(5*nr**2))
    f = 0
    k = 1
    do d = 1, 5
    do j = 1, nr
    do i = 1, nr
      f = f + 1
      mesh%boundary(2)%face(f)%element_id   = LinearElementIndex(d, i, j, k, nr)
      mesh%boundary(2)%face(f)%element_face = 5
    end do
    end do
    end do

    ! top
    mesh%boundary(3)%id = 3
    mesh%boundary(3)%name = 'top'
    allocate(mesh%boundary(3)%face(5*nr**2))
    f = 0
    k = nz
    do d = 1, 5
    do j = 1, nr
    do i = 1, nr
      f = f + 1
      mesh%boundary(3)%face(f)%element_id   = LinearElementIndex(d, i, j, k, nr)
      mesh%boundary(3)%face(f)%element_face = 6
    end do
    end do
    end do

    if (periodic) then
      mesh%boundary(2) % coupled  =  3
      mesh%boundary(2) % polarity = -1
      mesh%boundary(2) % map(3,4) =  h
      mesh%boundary(3) % coupled  =  2
      mesh%boundary(3) % polarity =  1
      mesh%boundary(3) % map(3,4) = -h
    end if

    ! vertex IDs ...............................................................

    call mesh % GenerateConsistentVertexIDs()

  end subroutine CreateCylinder

  !-----------------------------------------------------------------------------
  !> Linear index of vertex i,j,k in subdomain d

  pure integer function LinearVertexIndex(d, i, j, k, nr) result(l)
    integer, intent(in) :: d  !< subdomain, 1 <= d <= 5
    integer, intent(in) :: i  !< first index
    integer, intent(in) :: j  !< second index
    integer, intent(in) :: k  !< third index (axial or z direction)
    integer, intent(in) :: nr !< number of intervals along the inner quadrangle

    ! linear index within cross section ........................................

    select case(d)
    case(1) ! subdomain 1: inner square
      l = l1(i,j)
    case(2) ! subdomain 2: eastern ring section
      if (i == 0) then ! vertex belongs to inner square
        l = l1(nr,j)
      else if (j == 0) then ! vertex belongs to subdomain 5
        l = l5(i,nr)
      else
        l = l2(i,j)
      end if
    case(3) ! subdomain 3: northern ring section
      if (i == 0) then ! vertex belongs to inner square
        l = l1(nr-j,nr)
      else if (j == 0) then ! vertex belongs to subdomain 2
        l = l2(i,nr)
      else
        l = l3(i,j)
      end if
    case(4) ! subdomain 4: western ring section
      if (i == 0) then ! vertex belongs to inner square
        l = l1(0,nr-j)
      else if (j == 0) then ! vertex belongs to subdomain 3
        l = l3(i,nr)
      else
        l = l4(i,j)
      end if
    case(5) ! subdomain 5: southern ring section
      if (i == 0) then ! vertex belongs to inner square
        l = l1(j,0)
      else if (j == 0) then ! vertex belongs to subdomain 4
        l = l4(i,nr)
      else
        l = l5(i,j)
      end if
    case default
      l = 0
    end select

    ! add axial offset .........................................................

    l = l + ((nr+1)**2 + 4 * nr**2) * k

  contains

    ! subdomain 1:  0 <= i,j <= nr
    pure integer function l1(i, j)
      integer, intent(in) :: i, j
      l1 = i + (nr+1)*j + 1
    end function l1

    ! subdomain 2:  1 <= i,j <= nr
    pure integer function l2(i, j)
      integer, intent(in) :: i, j
      l2 = i + nr*(j-1) + (nr+1)**2
    end function l2

    ! subdomain 3:  1 <= i,j <= nr
    pure integer function l3(i, j)
      integer, intent(in) :: i, j
      l3 = i + nr*(j-1) + (nr+1)**2 + nr**2
    end function l3

    ! subdomain 4:  1 <= i,j <= nr
    pure integer function l4(i, j)
      integer, intent(in) :: i, j
      l4 = i + nr*(j-1) + (nr+1)**2 + 2 * nr**2
    end function l4

    ! subdomain 5:  1 <= i,j <= nr
    pure integer function l5(i, j)
      integer, intent(in) :: i, j
      l5 = i + nr*(j-1) + (nr+1)**2 + 3 * nr**2
    end function l5

  end function LinearVertexIndex

  !-----------------------------------------------------------------------------
  !> Linear index of element i,j,k in subdomain d

  pure integer function LinearElementIndex(d, i, j, k, nr) result(l)
    integer, intent(in) :: d  !< subdomain, 1 <= d <= 5
    integer, intent(in) :: i  !< first index
    integer, intent(in) :: j  !< second index
    integer, intent(in) :: k  !< third index (axial or z direction)
    integer, intent(in) :: nr  !< number of intervals along the inner quadrangle

    ! linear index within cross section ........................................

    select case(d)
    case(1) ! subdomain 1: inner square
      l = i + nr*(j-1)
    case(2) ! subdomain 2: eastern ring section
      l = i + nr*(j-1) + nr**2
    case(3) ! subdomain 3: northern ring section
      l = i + nr*(j-1) + 2 * nr**2
    case(4) ! subdomain 4: western ring section
      l = i + nr*(j-1) + 3 * nr**2
    case(5) ! subdomain 5: southern ring section
      l = i + nr*(j-1) + 4 * nr**2
    case default
      l = 0
    end select

    ! add axial offset .........................................................

    l = l + 5 * nr**2 * (k-1)

  end function LinearElementIndex

  !-----------------------------------------------------------------------------
  !> Point coordinates in the cylinder with radius 1 and axis (∓1,0,0)

  pure function PointCoordinates(d, xi, eta, zeta, r, h) result(x)
    integer,   intent(in) :: d    !< subdomain, 1 <= d <= 5
    real(RNP), intent(in) :: xi   !< first coordinate,  -1 <= xi   <= 1
    real(RNP), intent(in) :: eta  !< second coordinate, -1 <= eta  <= 1
    real(RNP), intent(in) :: zeta !< third coordinate,  -1 <= zeta <= 1
    real(RNP), intent(in) :: r    !< cylinder radius
    real(RNP), intent(in) :: h    !< cylinder height
    real(RNP) :: x(3)

    ! local variables ..........................................................

    real(RNP) :: a, b, c0, c1

    a = RI * sin(PI/12)
    b = ONE / sqrt(TWO)

    c0  = (b + a) * HALF
    c1  = (b - a) * HALF

    ! position in the square [-1,1]^2 ..........................................

    select case(d)

    case(1) ! subdomain 1: inner quadrangle
      x = ( AIW(-eta) * (1 -  xi)  +  AIE( eta) * (1 +  xi)   &
          + AIS(  xi) * (1 - eta)  +  AIN(- xi) * (1 + eta)   &
          ) / 2                                               &
        - ( [ -a, -a, zeta] * (1 - xi) * (1 - eta)            &
          + [  a, -a, zeta] * (1 + xi) * (1 - eta)            &
          + [ -a,  a, zeta] * (1 - xi) * (1 + eta)            &
          + [  a,  a, zeta] * (1 + xi) * (1 + eta)            &
          ) / 4

    case(2) ! subdomain 2: eastern ring section
      x = ( AIE( eta) * (1 -  xi)  +  AOE( eta) * (1 +  xi)   &
          + RSE(  xi) * (1 - eta)  +  RNE(  xi) * (1 + eta)   &
          ) / 2                                               &
        - ( [  a, -a, zeta] * (1 - xi) * (1 - eta)            &
          + [  b, -b, zeta] * (1 + xi) * (1 - eta)            &
          + [  a,  a, zeta] * (1 - xi) * (1 + eta)            &
          + [  b,  b, zeta] * (1 + xi) * (1 + eta)            &
          ) / 4

    case(3) ! subdomain 3: northern ring section
      x = ( AIN( eta) * (1 -  xi)  +  AON( eta) * (1 +  xi)   &
          + RNE(  xi) * (1 - eta)  +  RNW(  xi) * (1 + eta)   &
          ) / 2                                               &
        - ( [  a,  a, zeta] * (1 - xi) * (1 - eta)            &
          + [  b,  b, zeta] * (1 + xi) * (1 - eta)            &
          + [ -a,  a, zeta] * (1 - xi) * (1 + eta)            &
          + [ -b,  b, zeta] * (1 + xi) * (1 + eta)            &
          ) / 4

    case(4) ! subdomain 4: western ring section
      x = ( AIW( eta) * (1 -  xi)  +  AOW( eta) * (1 +  xi)   &
          + RNW(  xi) * (1 - eta)  +  RSW(  xi) * (1 + eta)   &
          ) / 2                                               &
        - ( [ -a,  a, zeta] * (1 - xi) * (1 - eta)            &
          + [ -b,  b, zeta] * (1 + xi) * (1 - eta)            &
          + [ -a, -a, zeta] * (1 - xi) * (1 + eta)            &
          + [ -b, -b, zeta] * (1 + xi) * (1 + eta)            &
          ) / 4

    case(5) ! subdomain 5: southern ring section
      x = ( AIS( eta) * (1 -  xi)  +  AOS( eta) * (1 +  xi)   &
          + RSW(  xi) * (1 - eta)  +  RSE(  xi) * (1 + eta)   &
          ) / 2                                               &
        - ( [ -a, -a, zeta] * (1 - xi) * (1 - eta)            &
          + [ -b, -b, zeta] * (1 + xi) * (1 - eta)            &
          + [  a, -a, zeta] * (1 - xi) * (1 + eta)            &
          + [  b, -b, zeta] * (1 + xi) * (1 + eta)            &
          ) / 4

    case default
      x = 0
    end select

    ! scale to cylinder radius and height
    x(1) = r   * x(1)
    x(2) = r   * x(2)
    x(3) = h/2 * x(3)

  contains

    ! inner east circular arc
    pure function AIE(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(RI*cos(u*PI/12) - DI, RNP), real(RI*sin(u*PI/12), RNP), zeta ]
    end function AIE

    ! inner north circular arc
    pure function AIN(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(-RI*sin(u*PI/12), RNP), real(RI*cos(u*PI/12) - DI, RNP), zeta ]
    end function AIN

    ! inner west circular arc
    pure function AIW(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(-RI*cos(u*PI/12) + DI, RNP), real(-RI*sin(u*PI/12), RNP), zeta ]
    end function AIW

    ! inner south circular arc
    pure function AIS(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(RI*sin(u*PI/12), RNP), real(-RI*cos(u*PI/12) + DI, RNP), zeta  ]
    end function AIS

    ! outer east circular arc
    pure function AOE(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(cos(u*PI/4), RNP), real(sin(u*PI/4), RNP), zeta  ]
    end function AOE

    ! outer north circular arc
    pure function AON(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(-sin(u*PI/4), RNP), real(cos(u*PI/4), RNP), zeta ]
    end function AON

    ! outer west circular arc
    pure function AOW(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(-cos(u*PI/4), RNP), real(-sin(u*PI/4), RNP), zeta ]
    end function AOW

    ! outer south circular arc
    pure function AOS(u) result(x)
      real(RNP), intent(in) :: u
      real(RNP) :: x(3)
      x = [ real(sin(u*PI/4), RNP), real(-cos(u*PI/4), RNP), zeta ]
    end function AOS

    ! north-east radius
    pure function RNE(v) result(x)
      real(RNP), intent(in) :: v
      real(RNP) :: x(3)
      x = [ (c0 + c1*v) * [ 1, 1],  zeta  ]
    end function RNE

    ! north-west radius
    pure function RNW(v) result(x)
      real(RNP), intent(in) :: v
      real(RNP) :: x(3)
      x = [ (c0 + c1*v) * [ -1, 1],  zeta  ]
    end function RNW

    ! south-west radius
    pure function RSW(v) result(x)
      real(RNP), intent(in) :: v
      real(RNP) :: x(3)
      x = [ (c0 + c1*v) * [ -1, -1],  zeta  ]
    end function RSW

    ! south-east radius
    pure function RSE(v) result(x)
      real(RNP), intent(in) :: v
      real(RNP) :: x(3)
      x = [ (c0 + c1*v) * [ 1, -1],  zeta  ]
    end function RSE

  end function PointCoordinates

  !=============================================================================

end submodule MP_CreateCylinder
