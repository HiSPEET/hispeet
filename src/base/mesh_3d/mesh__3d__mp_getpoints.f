!> summary:  Generation of mesh points to given order and basis type
!> author:   Joerg Stiller
!> date:     2020/11/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_GetPoints
  use Constants, only: HALF
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generates Gauss-Lobatto or Gauss points to all elements of a mesh partition
  !>
  !> The routine provides the element points for one of the following bases:
  !>
  !>   -  Lagrange polynomials to Gauss-Legendre points (nodes = 'G')
  !>   -  Lagrange polynomials to Gauss-Lobatto-Legendre points (nodes = 'L')

  module subroutine GetPoints(mesh, po, nodes, x)
    class(Mesh_3D),          intent(in)  :: mesh  !< mesh parition
    integer,                 intent(in)  :: po    !< polynomial order
    character,     optional, intent(in)  :: nodes !< 'G' or 'L' ['L']
    real(RNP),  allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points

    type(StandardElementOperators_1D) :: sop

    !$omp single
    allocate(x(0:po, 0:po, 0:po, mesh%n_elem, 3))
    !$omp end single

    if (mesh % n_elem < 1) return

    sop = StandardElementOperators_1D(po, nodes, no_vdm = .true.)

    if (mesh % regular) then
      call GetRegularMeshPoints(mesh, sop, x)
    else
      call GetDeformedMeshPoints(mesh, sop, x)
    end if

  end subroutine GetPoints

  !-----------------------------------------------------------------------------
  !> Generates points for a regular mesh

  subroutine GetRegularMeshPoints(mesh, sop, x)
    class(Mesh_3D),                     intent(in) :: mesh !< mesh parition
    class(StandardElementOperators_1D), intent(in) :: sop  !< standard operators
    real(RNP), contiguous, intent(out) :: x(0:,0:,0:,:,:)  !< mesh points

    real(RNP), dimension(0:sop%po)   :: x1, x2, x3
    real(RNP), dimension(1:sop%po-1) :: ys
    integer :: e, i, j, k

    ys = (sop % x(1:sop%po-1) + 1) / 2

    !$omp do
    do e = 1, mesh % n_elem
      associate(geo => mesh % element(e) % geometry)

        ! boundaries
        x1(      0 ) = geo % x_e(      0,      0,      0, 1)
        x1( sop%po ) = geo % x_e( geo%po,      0,      0, 1)
        x2(      0 ) = geo % x_e(      0,      0,      0, 2)
        x2( sop%po ) = geo % x_e(      0, geo%po,      0, 2)
        x3(      0 ) = geo % x_e(      0,      0,      0, 3)
        x3( sop%po ) = geo % x_e(      0,      0, geo%po, 3)

        ! 1D point distributions
        x1( 1:sop%po-1 ) = x1(0) + (x1(sop%po) - x1(0)) * ys
        x2( 1:sop%po-1 ) = x2(0) + (x2(sop%po) - x2(0)) * ys
        x3( 1:sop%po-1 ) = x3(0) + (x3(sop%po) - x3(0)) * ys

        ! element points
        do k = 0, sop%po
        do j = 0, sop%po
        do i = 0, sop%po
          x(i,j,k,e,1) = x1(i)
          x(i,j,k,e,2) = x2(j)
          x(i,j,k,e,3) = x3(k)
        end do
        end do
        end do

      end associate
    end do

  end subroutine GetRegularMeshPoints

  !-----------------------------------------------------------------------------
  !> Generates points for a nonuniform mesh

  subroutine GetDeformedMeshPoints(mesh, sop, x)
    class(Mesh_3D),                     intent(in) :: mesh !< mesh parition
    class(StandardElementOperators_1D), intent(in) :: sop  !< standard operators
    real(RNP), contiguous, intent(out) :: x(0:,0:,0:,:,:)  !< mesh points

    type(StandardElementOperators_1D)      :: gop(1:mesh%p_geom)
    type(EmbeddedInterpolationOperator_1D) :: iop(1:mesh%p_geom)
    integer :: d, e, pg, pm

    pm = sop % po

    !$omp do
    do e = 1, mesh % n_elem
      associate(geo => mesh % element(e) % geometry)

        pg = geo % po

        if (pm == pg) then
          do d = 1, 3
            x(:,:,:,e,d) = geo % x_e(:,:,:,d)
          end do

        else

          if (gop(pg) % po /= pg) then
            ! initialize interpolation operator
            gop(pg) = StandardElementOperators_1D(pg, nodes = 'L', no_vdm = .true.)
            iop(pg) = EmbeddedInterpolationOperator_1D(gop(pg), sop%x)
          end if

          do d = 1, 3
            call InterpolateCoords( pg, pm, iop(pg) % A     &
                                  , xg = geo % x_e(:,:,:,d) &
                                  , xm = x(:,:,:,e,d)       )
          end do

        end if

      end associate
    end do

  end subroutine GetDeformedMeshPoints

  !-----------------------------------------------------------------------------
  !> Interpolates element mesh point coordinates

  subroutine InterpolateCoords(pg, pm, A, xg, xm)
    integer,   intent(in)  :: pg
    integer,   intent(in)  :: pm
    real(RNP), intent(in)  :: A  (0:pm, 0:pg)
    real(RNP), intent(in)  :: xg (0:pg, 0:pg, 0:pg)
    real(RNP), intent(out) :: xm (0:pm, 0:pm, 0:pm)

    real(RNP) :: z1(0:pm, 0:pg, 0:pg)
    real(RNP) :: z2(0:pm, 0:pm, 0:pg)

    integer :: i, j, k, p

    do k = 0, pg
    do j = 0, pg
    do i = 0, pm
      z1(i,j,k) = 0
      do p = 0, pg
        z1(i,j,k) = z1(i,j,k) + A(i,p) * xg(p,j,k)
      end do
    end do
    end do
    end do

    do k = 0, pg
    do j = 0, pm
    do i = 0, pm
      z2(i,j,k) = 0
      do p = 0, pg
        z2(i,j,k) = z2(i,j,k) + A(j,p) * z1(i,p,k)
      end do
    end do
    end do
    end do

    do k = 0, pm
    do j = 0, pm
    do i = 0, pm
      xm(i,j,k) = 0
      do p = 0, pg
        xm(i,j,k) = xm(i,j,k) + A(k,p) * z2(i,j,p)
      end do
    end do
    end do
    end do

  end subroutine InterpolateCoords

  !=============================================================================

end submodule MP_GetPoints
