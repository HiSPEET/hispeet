!> summary:  Cartesian mesh partition
!> author:   Joerg Stiller
!> date:     2017/04/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh partition
!===============================================================================

submodule(CART__Mesh_Partition) CART__Mesh_Partition__Get_Points
  use Kind_Parameters, only: RNP
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
!>
!> GLL is the default, but po = 0 implies GL irrespective of the chosen basis.

module subroutine GetPoints(mesh, po, x, basis)

  class(MeshPartition),   intent(in)  :: mesh         !< mesh parition
  integer,                intent(in)  :: po           !< polynomial order
  real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points
  character(len=*),       intent(in)  :: basis        !< set: 'GL' or 'GLL'
  optional :: basis

  type(StandardOperators1D) :: standard_op

  real(RNP), allocatable :: x1(:), x2(:), x3(:)
  real(RNP) :: dx(3)
  integer   :: e, i, j, k

  ! initialization .............................................................

  allocate(x(0:po, 0:po, 0:po, mesh%ne, 3))
  allocate(x1(0:po), x2(0:po), x3(0:po))

  call standard_op % New(po, basis)

  ! create element points ......................................................

  associate(xs => standard_op % x)

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

end submodule CART__Mesh_Partition__Get_Points
