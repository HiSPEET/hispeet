!> summary:  Maps shell collocation points from GMSH to Cartesian array
!> author:   Matthias Frey
!> date:     2023/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Import_GMSH__3D) MP_MapShellNodes
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Maps shell vertex, edge and surface nodes to collocation points

  module subroutine MapShellNodes(v, e, s, shellOrder, points)
    integer, intent(in)  :: shellOrder !< shell order
    integer, intent(in)  :: v(:)       !< vertex nodes
    integer, intent(in)  :: e(:)       !< edge nodes without end nodes
    integer, intent(in)  :: s(:)       !< surface nodes without boundary
    integer, allocatable, intent(out) :: points(:,:,:) !< collocation points

    ! auxiliary variables ......................................................

    integer :: i, j, k, l, m, n, r
    integer :: np, np_edge, np_face
    integer :: numShells, numCircles, circleOrder

    integer, allocatable :: p_face(:,:)
    integer, allocatable :: p_circle(:,:)
    integer, allocatable :: p_edge(:,:,:,:)

    ! prerequisites ............................................................

    np         = shellOrder + 1           ! number of nodes per direction
    np_edge    = shellOrder - 1           ! number of nodes per edge
    np_face    = np_edge**2               ! number of nodes per face
    numShells  = floor(shellOrder/2.)     ! number of shells
    numCircles = floor(shellOrder/2.) - 1 ! number of circles

    allocate(points(np,np,np), source = 0)
    allocate(p_face(6,np_face), source = 0)
    allocate(p_edge(6,numShells,4,shellOrder-1), source = 0)

    ! vertex points, assuming rotational numbering .............................

    points(np, 1, 1) = v(1)
    points(np,np, 1) = v(2)
    points(np,np,np) = v(3)
    points(np, 1,np) = v(4)
    points( 1, 1, 1) = v(5)
    points( 1,np, 1) = v(6)
    points( 1,np,np) = v(7)
    points( 1, 1,np) = v(8)

    ! edge points, assuming rotational numbering ...............................

    points(np,2:np-1,1)     = e( 0*np_edge+1: 1*np_edge) ! edge  1 (vertex 1->2)
    points(np,1,2:np-1)     = e( 1*np_edge+1: 2*np_edge) ! edge  2 (vertex 1->4)
    points(np-1:2:-1,1,1)   = e( 2*np_edge+1: 3*np_edge) ! edge  3 (vertex 1->5)
    points(np,np,2:np-1)    = e( 3*np_edge+1: 4*np_edge) ! edge  4 (vertex 2->3)
    points(np-1:2:-1,np,1)  = e( 4*np_edge+1: 5*np_edge) ! edge  5 (vertex 2->6)
    points(np,np-1:2:-1,np) = e( 5*np_edge+1: 6*np_edge) ! edge  6 (vertex 3->4)
    points(np-1:2:-1,np,np) = e( 6*np_edge+1: 7*np_edge) ! edge  7 (vertex 3->7)
    points(np-1:2:-1,1,np)  = e( 7*np_edge+1: 8*np_edge) ! edge  8 (vertex 4->8)
    points(1,2:np-1,1)      = e( 8*np_edge+1: 9*np_edge) ! edge  9 (vertex 5->6)
    points(1,1,2:np-1)      = e( 9*np_edge+1:10*np_edge) ! edge 10 (vertex 5->8)
    points(1,np,2:np-1)     = e(10*np_edge+1:11*np_edge) ! edge 11 (vertex 6->7)
    points(1,np-1:2:-1,np)  = e(11*np_edge+1:12*np_edge) ! edge 12 (vertex 7->8)

    ! surface points, assuming rotational numbering ............................

    do l = 1, 6

      ! map surface nodes to six faces
      p_face(l,:) = s((l-1)*np_face+1:l*np_face)

      ! map surface nodes to circle nodes
      circleOrder = shellOrder - 2
      k = 1
      do i = 1, numCircles
        allocate(p_circle(6,np_circle(circleOrder)), source = 0)

        p_circle(l,:) = p_face(l,k:k+np_circle(circleOrder)-1)

        m = 2
        n = 5
        do j = 1,4
          if(j == 4) m = 1

          ! map circle nodes to face edges
          p_edge(l,i,j,1:circleOrder+1) = &
                   [p_circle(l,j), p_circle(l,n:n+circleOrder-2), p_circle(l,m)]

          m = m + 1
          n = n + circleOrder-2 + 1
        end do

        k = k + np_circle(circleOrder)
        circleOrder = circleOrder - 2
        deallocate(p_circle)
      end do

      ! map center of shell faces
      select case(mod(shellOrder,2))
      case(0)
        p_edge(l,numShells,1:4,:) = p_face(l,np_face)

      case(1)
        m = size(p_face(l,:))
        p_edge(l,numShells,1,1:2) = [ p_face(l,m-3), p_face(l,m-2) ]
        p_edge(l,numShells,2,1:2) = [ p_face(l,m-2), p_face(l,m-1) ]
        p_edge(l,numShells,3,1:2) = [ p_face(l,m-1), p_face(l,m  ) ]
        p_edge(l,numShells,4,1:2) = [ p_face(l,m  ), p_face(l,m-3) ]
      end select

    end do

    ! map face circles to each surface recursively .............................

    ! surface 1 (vertex 1,2,3,4 rotational)
    n = np_edge
    do r = 1, numCircles + 1
      points(np,1+r,1+r:np-r)     = p_edge(1,r,1,1:n)
      points(np,1+r:np-r,np-r)    = p_edge(1,r,2,1:n)
      points(np,np-r,np-r:1+r:-1) = p_edge(1,r,3,1:n)
      points(np,np-r:1+r:-1,1+r)  = p_edge(1,r,4,1:n)
      n = n-2
    end do

    ! surface 2 (vertex 1,2,3,4 rotational)
    n = np_edge
    do r = 1, numCircles + 1
      points(np-r,1+r:np-r,1)    = p_edge(2,r,1,1:n)
      points(np-r:1+r:-1,np-r,1) = p_edge(2,r,2,1:n)
      points(1+r,np-r:1+r:-1,1)  = p_edge(2,r,3,1:n)
      points(1+r:np-r,1+r,1)     = p_edge(2,r,4,1:n)
      n = n-2
    end do

    ! surface 3 (vertex 1,5,8,4 rotational)
    n = np_edge
    do r =1, numCircles + 1
      points(np-r:1+r:-1,1,1+r)  = p_edge(3,r,1,1:n)
      points(1+r,1,1+r:np-r)     = p_edge(3,r,2,1:n)
      points(1+r:np-r,1,np-r)    = p_edge(3,r,3,1:n)
      points(np-r,1,np-r:1+r:-1) = p_edge(3,r,4,1:n)
      n = n-2
    end do

    ! surface 4 (vertex 2,3,7,6 rotational)
    n = np_edge
    do r = 1, numCircles + 1
      points(np-r,np,1+r:np-r)    = p_edge(4,r,1,1:n)
      points(np-r:1+r:-1,np,np-r) = p_edge(4,r,2,1:n)
      points(1+r,np,np-r:1+r:-1)  = p_edge(4,r,3,1:n)
      points(1+r:np-r,np,1+r)     = p_edge(4,r,4,1:n)
      n = n-2
    end do

    ! surface 5 (vertex 4,3,7,8 rotational)
    n = np_edge
    do r = 1, numCircles + 1
      points(np-r,np-r:1+r:-1,np) = p_edge(5,r,1,1:n)
      points(np-r:1+r:-1,1+r,np)  = p_edge(5,r,2,1:n)
      points(1+r,1+r:np-r,np)     = p_edge(5,r,3,1:n)
      points(1+r:np-r,np-r,np)    = p_edge(5,r,4,1:n)
      n = n-2
    end do

    ! surface 6 (vertex 5,6,7,8 rotational)
    n = np_edge
    do r = 1, numCircles + 1
      points(1,1+r:np-r,1+r)     = p_edge(6,r,1,1:n)
      points(1,np-r,1+r:np-r)    = p_edge(6,r,2,1:n)
      points(1,np-r:1+r:-1,np-r) = p_edge(6,r,3,1:n)
      points(1,1+r,np-r:1+r:-1)  = p_edge(6,r,4,1:n)
      n = n-2
    end do

  end subroutine MapShellNodes

  !-----------------------------------------------------------------------------
  !> Returns number of points per circle for a given circle order

  integer function np_circle(circleOrder)
    integer, intent(in) :: circleOrder !< circle order

    np_circle = (circleOrder+1)**2 - (circleOrder-1)**2

  end function

  !=============================================================================

end submodule MP_MapShellNodes
