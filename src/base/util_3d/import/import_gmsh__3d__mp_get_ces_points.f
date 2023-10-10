!> summary:  Returns a 3D-matrix of hexa element points of a GMSH .msh file
!> author:   Matthias Frey
!> date:     2023/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Import_GMSH__3D) MP_Get_CES_Points
  implicit none

contains
      
  !-----------------------------------------------------------------------------
  !> Takes corner, edge and surfaces points of any mesh order and 
  !> returns a 3D-matrix with corresponding node tags of cube layer
      

  module subroutine getCESpoints(v,e,s,shellOrder,p)
    !implicit none
    integer, intent(in)  :: shellOrder                !< shell order
    integer, intent(in)  :: v(:), e(:), s(:)          !< points of vertices, edges, surfaces
    integer, allocatable, intent(out) :: p(:,:,:)     !< 3D point array to be returned

    ! auxiliary variables ......................................................

    integer :: i, j, k, l, m, n, r
    integer :: np_dir, np_edge, np_surface
    integer :: numShells, numCircles, circleOrder

    integer, allocatable :: p_surf(:,:)        !< array to subdivide s vector into six sufaces    
    integer, allocatable :: p_circle(:,:)      !< array to subdivide each p_surf vector into circles
    integer, allocatable :: p_cl_edge(:,:,:,:) !< surface circle edge (number, circle, circular row, points)

  

    ! prerequisites ............................................................
  
    np_dir     = shellOrder+1                ! number of points per direction (edge with vertices included)
    np_edge    = shellOrder-1                ! number of points per edge (without corner vertices)
    np_surface = np_edge**2                  ! number of points per surface
    numShells  = floor(shellOrder/2.)        ! number of shells
    numCircles = floor(shellOrder/2.) - 1    ! number of circles

    allocate(p(np_dir,np_dir,np_dir), source = 0)
    allocate(p_surf(6,np_surface), source = 0)
    allocate(p_cl_edge(6,numShells,4,shellOrder-1), source = 0)


    ! allocate 8 corner vertices (rotational numbering) ........................

    p(np_dir,1,1)           = v(1) 
    p(np_dir,np_dir,1)      = v(2)
    p(np_dir,np_dir,np_dir) = v(3) 
    p(np_dir,1,np_dir)      = v(4) 
    p(1,1,1)                = v(5)  
    p(1,np_dir,1)           = v(6) 
    p(1,np_dir,np_dir)      = v(7) 
    p(1,1,np_dir)           = v(8)


    ! allocate edge points (rotational numbering) ..............................

    p(np_dir,2:np_dir-1,1)         = e( 0*np_edge+1: 1*np_edge)  ! edge I    (vertex 1->2)
    p(np_dir,1,2:np_dir-1)         = e( 1*np_edge+1: 2*np_edge)  ! edge II   (vertex 1->4)
    p(np_dir-1:2:-1,1,1)           = e( 2*np_edge+1: 3*np_edge)  ! edge III  (vertex 1->5)
    p(np_dir,np_dir,2:np_dir-1)    = e( 3*np_edge+1: 4*np_edge)  ! edge IV   (vertex 2->3)
    p(np_dir-1:2:-1,np_dir,1)      = e( 4*np_edge+1: 5*np_edge)  ! edge V    (vertex 2->6)
    p(np_dir,np_dir-1:2:-1,np_dir) = e( 5*np_edge+1: 6*np_edge)  ! edge VI   (vertex 3->4)
    p(np_dir-1:2:-1,np_dir,np_dir) = e( 6*np_edge+1: 7*np_edge)  ! edge VII  (vertex 3->7)
    p(np_dir-1:2:-1,1,np_dir)      = e( 7*np_edge+1: 8*np_edge)  ! edge VIII (vertex 4->8)
    p(1,2:np_dir-1,1)              = e( 8*np_edge+1: 9*np_edge)  ! edge IX   (vertex 5->6)
    p(1,1,2:np_dir-1)              = e( 9*np_edge+1:10*np_edge)  ! edge X    (vertex 5->8)
    p(1,np_dir,2:np_dir-1)         = e(10*np_edge+1:11*np_edge)  ! edge XI   (vertex 6->7)
    p(1,np_dir-1:2:-1,np_dir)      = e(11*np_edge+1:12*np_edge)  ! edge XII  (vertex 7->8)


    ! allocate surface points ..................................................

    do l = 1,6
      
      ! subdivide and map the surface point vector s into six surfaces p_surf
      p_surf(l,:) = s((l-1)*np_surface+1:l*np_surface)

      ! subdivide and map each surface vector p_surf into ncircles circles p_circle
      circleOrder = shellOrder - 2
      k = 1
      do i = 1, numCircles
        allocate(p_circle(6,np_circle(circleOrder)), source = 0) 

        p_circle(l,:) = p_surf(l,k:k+np_circle(circleOrder)-1)

        m = 2
        n = 5
        do j = 1,4
          if(j == 4) m = 1

          ! subdivide and map each circle vector to 4 surface edges p_cl_edge
          p_cl_edge(l,i,j,1:circleOrder+1) = [p_circle(l,j), p_circle(l,n:n+circleOrder-2), p_circle(l,m)]
       
          m = m+1
          n = n + circleOrder-2 + 1
        end do

        k = k + np_circle(circleOrder)
        circleOrder = circleOrder - 2
        deallocate(p_circle) 
      end do

      ! map inner incomplete round of each surface
      select case(mod(shellOrder,2))
      case(0)
      p_cl_edge(l,numShells,1:4,:) = p_surf(l,size(p_surf(l,:)))

      case(1)
      p_cl_edge(l,numShells,1,:)   =  [p_surf(l,size(p_surf(l,:))-3), p_surf(l,size(p_surf(l,:))-2)]
      p_cl_edge(l,numShells,2,:)   =  [p_surf(l,size(p_surf(l,:))-2), p_surf(l,size(p_surf(l,:))-1)]
      p_cl_edge(l,numShells,3,:)   =  [p_surf(l,size(p_surf(l,:))-1), p_surf(l,size(p_surf(l,:)))]
      p_cl_edge(l,numShells,4,:)   =  [p_surf(l,size(p_surf(l,:))),   p_surf(l,size(p_surf(l,:))-3)]
      end select

    end do




    ! map surface points to each surface recursively ..........................
   
    ! surface I (vertex 1,2,3,4 rotational) 
    n = np_edge
    do r = 1, numShells
      p(np_dir,1+r,1+r:np_dir-r)          = p_cl_edge(1,r,1,1:n)
      p(np_dir,1+r:np_dir-r,np_dir-r)     = p_cl_edge(1,r,2,1:n)
      p(np_dir,np_dir-r,np_dir-r:1+r:-1)  = p_cl_edge(1,r,3,1:n)
      p(np_dir,np_dir-r:1+r:-1,1+r)       = p_cl_edge(1,r,4,1:n)
      n = n-2
    end do

    ! surface II (vertex 1,2,3,4 rotational) 
    n = np_edge
    do r = 1, numShells
      p(np_dir-r,1+r:np_dir-r,1)          = p_cl_edge(2,r,1,1:n)
      p(np_dir-r:1+r:-1,np_dir-r,1)       = p_cl_edge(2,r,2,1:n)
      p(1+r,np_dir-r:1+r:-1,1)            = p_cl_edge(2,r,3,1:n)
      p(1+r:np_dir-r,1+r,1)               = p_cl_edge(2,r,4,1:n)
      n = n-2
    end do

    ! surface III (vertex 1,5,8,4 rotational) 
    n = np_edge
    do r =1, numShells
      p(np_dir-r:1+r:-1,1,1+r)            = p_cl_edge(3,r,1,1:n)
      p(1+r,1,1+r:np_dir-r)               = p_cl_edge(3,r,2,1:n)
      p(1+r:np_dir-r,1,np_dir-r)          = p_cl_edge(3,r,3,1:n)
      p(np_dir-r,1,np_dir-r:1+r:-1)       = p_cl_edge(3,r,4,1:n)
      n = n-2
    end do

    ! surface IV (vertex 2,3,7,6 rotational) 
    n = np_edge
    do r = 1, numShells
      p(np_dir-r,np_dir,1+r:np_dir-r)     = p_cl_edge(4,r,1,1:n)
      p(np_dir-r:1+r:-1,np_dir,np_dir-r)  = p_cl_edge(4,r,2,1:n)
      p(1+r,np_dir,np_dir-r:1+r:-1)       = p_cl_edge(4,r,3,1:n)
      p(1+r:np_dir-r,np_dir,1+r)          = p_cl_edge(4,r,4,1:n)
      n = n-2
    end do

    ! surface V (vertex 4,3,7,8 rotational) 
    n = np_edge
    do r = 1, numShells
      p(np_dir-r,np_dir-r:1+r:-1,np_dir)  = p_cl_edge(5,r,1,1:n)
      p(np_dir-r:1+r:-1,1+r,np_dir)       = p_cl_edge(5,r,2,1:n)
      p(1+r,1+r:np_dir-r,np_dir)          = p_cl_edge(5,r,3,1:n)
      p(1+r:np_dir-r,np_dir-r,np_dir)     = p_cl_edge(5,r,4,1:n)
      n = n-2
    end do

    ! surface VI (vertex 5,6,7,8 rotational) 
    n = np_edge
    do r = 1, numShells
      p(1,1+r:np_dir-r,1+r)               = p_cl_edge(6,r,1,1:n)
      p(1,np_dir-r,1+r:np_dir-r)          = p_cl_edge(6,r,2,1:n)
      p(1,np_dir-r:1+r:-1,np_dir-r)       = p_cl_edge(6,r,3,1:n)
      p(1,1+r,np_dir-r:1+r:-1)            = p_cl_edge(6,r,4,1:n)
      n = n-2
    end do


  end subroutine getCESpoints



  !> Function np_circle returns number of points per circle for a given circle order
    
  function np_circle(circleOrder) result(npoints)
  integer, intent (in) :: circleOrder          !< shell order
  integer              :: npoints              !< number of points per shell

  npoints = (circleOrder+1)**2 - (circleOrder-1)**2 
  end function

end submodule MP_Get_CES_Points
