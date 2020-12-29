!> summary:  Generation of mesh points to given order an basis type
!> author:   Joerg Stiller
!> date:     2020/11/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_GetPoints
  use Constants, only: HALF
  use Standard_Operators_1D
  use Embedded_Interpolation_3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generates Gauss-Lobatto or Gauss points to all elements of a mesh partition
  !>
  !> The routine provides the element points for one of the following bases:
  !>
  !>   *  Lagrange polynomials to Gauss-Legendre points (basis = 'GL')
  !>   *  Lagrange polynomials to Gauss-Lobatto-Legendre points (basis = 'GLL')

  module subroutine GetPoints(mesh, po, basis, x)
    class(Mesh3d_Partition), intent(in)  :: mesh         !< mesh parition
    integer,                 intent(in)  :: po           !< polynomial order
    character(len=*),        intent(in)  :: basis        !< 'GL' or 'GLL'
    real(RNP), allocatable,  intent(out) :: x(:,:,:,:,:) !< mesh points

    type(StandardOperators1D) :: sop

    sop = StandardOperators1D(po, basis, no_vdm = .true.)

    !$omp single
    allocate(x(0:po, 0:po, 0:po, mesh%n_elem, 3))
    !$omp end single

    if (mesh % n_elem < 1) return

    if (mesh % regular) then
      call GetRegularMeshPoints(mesh, sop, x)
    else
      call GetDeformedMeshPoints(mesh, sop, x)
    end if

  end subroutine GetPoints

  !-----------------------------------------------------------------------------
  !> Generates points for a regular mesh

  subroutine GetRegularMeshPoints(mesh, sop, x)
    class(Mesh3d_Partition),    intent(in)  :: mesh !< mesh parition
    class(StandardOperators1D), intent(in)  :: sop  !< standard operators
    real(RNP), contiguous,      intent(out) :: x(0:,0:,0:,:,:) !< mesh points

    real(RNP), dimension(0:sop%po)   :: x1, x2, x3
    real(RNP), dimension(1:sop%po-1) :: ys
    integer   :: e, i, j, k

    associate( pg => mesh % p_geom  &
             , po => sop  % po      &
             , xs => sop  % x       )

      ys = (sop % x(1:po-1) + 1) / 2

      !$omp do
      do e = 1, mesh % n_elem

        ! boundaries
        x1( 0) = mesh % x_elem( 0,  0,  0, e, 1)
        x1(po) = mesh % x_elem(po,  0,  0, e, 1)
        x2( 0) = mesh % x_elem( 0,  0,  0, e, 2)
        x2(po) = mesh % x_elem( 0, po,  0, e, 2)
        x3( 0) = mesh % x_elem( 0,  0,  0, e, 3)
        x3(po) = mesh % x_elem( 0,  0, po, e, 3)

        ! 1D point distributions
        x1(1:po-1) = x1(0) + (x1(po) - x1(0)) * ys
        x2(1:po-1) = x2(0) + (x2(po) - x2(0)) * ys
        x3(1:po-1) = x3(0) + (x3(po) - x3(0)) * ys

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

      end do

    end associate

  end subroutine GetRegularMeshPoints

  !-----------------------------------------------------------------------------
  !> Generates points for a nonuniform mesh

  subroutine GetDeformedMeshPoints(mesh, sop, x)
    class(Mesh3d_Partition),    intent(in)  :: mesh !< mesh parition
    class(StandardOperators1D), intent(in)  :: sop  !< standard operators
    real(RNP), contiguous,      intent(out) :: x(0:,0:,0:,:,:) !< mesh points

    type(StandardOperators1D) :: gop
    type(InterpolationOperator3D) :: iop

    gop = StandardOperators1D(mesh % p_geom, basis = 'GLL', no_vdm = .true.)
    iop = InterpolationOperator3D(gop, sop%x)

    call iop % Apply(mesh % x_elem, x)

  end subroutine GetDeformedMeshPoints

  !=============================================================================

end submodule MP_GetPoints
