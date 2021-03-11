!> summary:  Creation of a 3d generic mesh in an annular gap
!> author:   Susanne Stimpert, Joerg Stiller
!> date:     2015/05/08, revised 2021/01/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!==============================================================================

submodule(Generic_Mesh_3d) MP_CreateAnnularGap
  use Constants
  use Gauss_Jacobi
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Creates a generic mesh for an annular gap
  !>
  !> This routine generates an annular gap with inner radius `r0`, outer radius
  !> `r1` and height `h`. The mesh consists of `nr` elements in radial direction,
  !> `np` elements in azimuthal direction and `nz` elements in axial direction.
  !> The polynomial degree of the GLL points is `po`. If `periodic == .true.`,
  !> axial periodicity is set.
  !>
  !> The workflow to set the vertices and the GLL points is as following:
  !>   - First the (nz+1)*(nr+1)*np vertices are set equidistantly in the unit
  !>     cube [0, 1]^3. They are needed to set the GLL points in the subcubes.
  !>   - In every subcube (po+1)^3 GLL points are initialized. In the same loop
  !>     they are mapped to the cartesian coordinates in the gap.
  !>   - At the end the vertex coordinates are transformed to the gap.
  !>                                         __________
  !>                                        /          \
  !>                __________             /            \
  !>             nz/         /|           |  __________  |
  !>              /________ / |         nz| /    ___   \ |
  !>             |         |  |    --->   |/    /   \   \|
  !>           nr|         |  |           |-nr-(     )   |
  !>             |         | /             \    \___/   /  <\
  !>             |_________|/               \__________/    / np
  !>                  np                                 __/
  !>

  module subroutine CreateAnnularGap(mesh, r0, r1, h, nr, np, nz, po, periodic)
    class(GenericMesh3d), intent(out) :: mesh  !< cylindrical 3d mesh
    real(RNP), intent(in) :: r0  !< inner radius
    real(RNP), intent(in) :: r1  !< outer radius
    real(RNP), intent(in) :: h   !< height
    integer,   intent(in) :: nr  !< num elements in radial    (r)   direction
    integer,   intent(in) :: np  !< num elements in azimuthal (phi) direction
    integer,   intent(in) :: nz  !< num elements in axial     (z)   direction
    integer,   intent(in) :: po  !< polynomial order of mesh elements
    logical,   intent(in) :: periodic !< switch for axial periodicity

    ! local variables
    integer   :: nv, ne, nb       ! number of vertices, elements and boundaries
    integer   :: i, j, k, p, q, r ! loop indices
    integer   :: l                ! linear element index
    integer   :: v                ! linear vertex index
    integer   :: f                ! face counter
    integer   :: b                ! boundary ID

    real(RNP) :: xc(0:po)         ! GLL points in [0,1]
    real(RNP) :: x(3), x0(3)      ! position vectors
    real(RNP) :: cr, cp, cz       ! 1/np, 1/nr, 1/nz

    ! prerequisites ...........................................................

    ! number of vertices, elements and boundaries
    nv = (nr+1) * np * (nz+1) ! num vertices
    ne = nr * np * nz         ! num elements
    nb = 4                    ! num boundaries (bottom, inner, outer, top)

    cr = ONE / nr
    cp = ONE / np
    cz = ONE / nz

    ! GLL points in [0,1]
    xc = HALF * (ONE + GLL_Points(po))

    ! allocate mesh components
    allocate(mesh%vertex(nv), mesh%element(ne), mesh%boundary(nb))

    ! set lexical numbering
    mesh%numbering = LEXICAL_NUMBERING

    ! vertices in cartesian coordinates (unit cube) ............................

    do k = 0, nz
    do j = 0, np-1
    do i = 0, nr
      v = VertexIndex(i, j, k, nz, np)
      mesh%vertex(v)%x  = [ i*cr, j*cp, k*cz ]
    end do
    end do
    end do

    ! elements ...................................................................

    mesh%element%typ   = HEXAHEDRAL_ELEMENT
    mesh%element%basis = GAUSS_LOBATTO_BASIS
    mesh%element%order = po

    do k = 1, nz
    do j = 1, np
    do i = 1, nr
      l = ElementIndex(i, j, k, nr, np)

      mesh%element(l)%id = l

      ! vertices
      mesh%element(l)%vertex = [ VertexIndex(i-1, j-1, k-1, nr, np), &
                                 VertexIndex(i  , j-1, k-1, nr, np), &
                                 VertexIndex(i-1, j  , k-1, nr, np), &
                                 VertexIndex(i  , j  , k-1, nr, np), &
                                 VertexIndex(i-1, j-1, k  , nr, np), &
                                 VertexIndex(i  , j-1, k  , nr, np), &
                                 VertexIndex(i-1, j  , k  , nr, np), &
                                 VertexIndex(i  , j  , k  , nr, np)  ]

      ! element collocation points
      allocate(mesh%element(l)%x((po+1)**3, 3))
      x0 = mesh%vertex(mesh%element(l)%vertex(1)) % x ! coordinates of vertex 1
      do r = 0, po
      do p = 0, po
      do q = 0, po
        v = 1 + p + (po+1) * (q + (po+1) * r)     ! linear point index
        x = x0 + [ cr*xc(p), cp*xc(q), cz*xc(r) ] ! point in unit cube
        call MapToGap(x, r0, r1, h)               ! map to annular gap
        mesh%element(l) % x(v,1:3) = x
      end do
      end do
      end do

    end do
    end do
    end do

    ! map vertices from unit cube to annular gap ...............................

    do v = 1, nv
      call MapToGap(mesh%vertex(v)%x, r0, r1, h)
    end do

    ! boundaries ...............................................................

    b = 1
    mesh%boundary(b)%id = 2
    mesh%boundary(b)%name = 'inner'
    allocate(mesh%boundary(b)%face(np*nz))
    f = 0
    do k = 1, nz
    do j = 1, np
    do i = 1, 1
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(i, j, k, nr, np)
      mesh%boundary(b)%face(f)%element_face = 1
    end do
    end do
    end do

    b = 2
    mesh%boundary(b)%id = 3
    mesh%boundary(b)%name = 'outer'
    allocate(mesh%boundary(b)%face(np*nz))
    f = 0
    do k = 1, nz
    do j = 1, np
    do i = nr, nr
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(i, j, k, nr, np)
      mesh%boundary(b)%face(f)%element_face = 2
    end do
    end do
    end do

    b = 3
    mesh%boundary(b)%id = 1
    mesh%boundary(b)%name = 'bottom'
    allocate(mesh%boundary(b)%face(nr*np))
    f = 0
    do k = 1, 1
    do j = 1, np
    do i = 1, nr
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(i, j, k, nr, np)
      mesh%boundary(b)%face(f)%element_face = 5
    end do
    end do
    end do

    b = 4
    mesh%boundary(b)%id = 4
    mesh%boundary(b)%name = 'top'
    allocate(mesh%boundary(b)%face(nr*np))
    f = 0
    do k = nz, nz
    do j =  1, np
    do i =  1, nr
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(i, j, k, nr, np)
      mesh%boundary(b)%face(f)%element_face = 6
    end do
    end do
    end do

    ! in the periodic case identify the bottom and the top
    if (periodic) then
      mesh%boundary(3)%coupled = 4
      mesh%boundary(4)%coupled = 3
    end if

    ! vertex IDs ...............................................................

    call mesh % GenerateConsistentVertexIDs()

  end subroutine CreateAnnularGap

  !------------------------------------------------------------------------------
  !> Linear vertex index based on lexical ordering

  function VertexIndex(i, j, k, nr, np) result(v)
    integer, intent(in)  :: i, j, k  !< triple vertex index starting at [0,0,0]
    integer, intent(in)  :: nr       !< num elements in radial direction
    integer, intent(in)  :: np       !< num elements in azimuthal direction
    integer              :: v        !< linear vertex index

    ! local variable
    integer :: j1

    j1 = modulo(j, np) ! closing the ring
    v  = 1 + i + (nr+1) * (j1 + np * k)

  end function VertexIndex

  !-----------------------------------------------------------------------------
  !> Linear element index based on lexical ordering

  pure function ElementIndex(i, j, k, nr, np) result(l)
    integer, intent(in)  :: i, j, k  !< triple index element starting at [1,1,1]
    integer, intent(in)  :: nr       !< num elements in radial direction
    integer, intent(in)  :: np       !< num elements in azimuthal direction
    integer              :: l        !< linear element index

    l = i + nr * (j-1 + np * (k-1))

  end function ElementIndex

  !-----------------------------------------------------------------------------
  !>  Maps x from unit cube to annular gap
  !>
  !> This routine maps a point of a cube [0, 1]^3 to the fitting position in a
  !> cylindrical shell with thickness (r1 - r0). The x and y coordinates are
  !> transformed to polar coordinates and scaled in radial direction. The z
  !> direction is scaled to the height h.

  subroutine MapToGap(x, r0, r1, h)
    real(RNP), intent(inout) :: x(3)    !< point to be mapped
    real(RNP), intent(in)    :: r0, r1  !< inner and outer radius of the shell
    real(RNP), intent(in)    :: h       !< height of the shell

    ! local variable
    real(RNP) :: r, phi

    r   = r0 + (r1 - r0)*(x(1))
    phi = 2 * PI * x(2)
    x   = [ r * cos(phi), r * sin(phi), h * x(3) ]

  end subroutine MapToGap

  !=============================================================================

end submodule MP_CreateAnnularGap
