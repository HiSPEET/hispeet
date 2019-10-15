!> summary:  Generation of mesh points to given order an basis type
!> author:   Joerg Stiller
!> date:     2017/04/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generation of mesh points to given order an basis type
!===============================================================================

submodule(CART__Mesh_Partition) MP_GetPoints
!  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Standard_Operators_1D
  implicit none

contains

!-------------------------------------------------------------------------------
!> Generates Gauss-Lobatto or Gauss points to all elements of a mesh partition
!>
!> The routine provides the element points for one of the following bases:
!>
!>   *  Lagrange polynomials to Gauss-Legendre points (basis = 'GL')
!>   *  Lagrange polynomials to Gauss-Lobatto-Legendre points (basis = 'GLL')

module subroutine GetPoints(mesh, po, basis, x)
  class(MeshPartition),   intent(in)  :: mesh         !< mesh parition
  integer,                intent(in)  :: po           !< polynomial order
  character(len=*),       intent(in)  :: basis        !< 'GL' or 'GLL'
  real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points

  type(StandardOperators1D) :: sop

  real(RNP), allocatable :: x1(:), x2(:), x3(:)
  real(RNP) :: dx(3)
  integer   :: e, i, j, k

  ! initialization .............................................................

  allocate(x(0:po, 0:po, 0:po, mesh%ne, 3))
  allocate(x1(0:po), x2(0:po), x3(0:po))

  sop = StandardOperators1D(po, basis, no_vdm = .true.)

  ! create element points ......................................................

  associate(xs => sop % x)

    do e = 1, mesh%ne

      associate(vertex => mesh%element(e)%vertex)

        dx = vertex(8)%x - vertex(1)%x

        ! 1D point distributions
        x1 = vertex(1)%x(1) + HALF * dx(1) * (xs + 1)
        x2 = vertex(1)%x(2) + HALF * dx(2) * (xs + 1)
        x3 = vertex(1)%x(3) + HALF * dx(3) * (xs + 1)

        ! element points
        do k = 0, po
        do j = 0, po
        do i = 0, po
          x(i,j,k,e,1) = x1(i)
          x(i,j,k,e,2) = x2(j)
          x(i,j,k,e,3) = x3(k)
        end do
        end do
        end do

      end associate

    end do

  end associate

end subroutine GetPoints

!===============================================================================

end submodule MP_GetPoints
