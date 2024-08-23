!> summary:  Generation approximate cuboids to mesh elements
!> author:   Joerg Stiller
!> date:     2022/06/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildCuboids
  use Standard_Element_Operators__1D
  use Element_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of approximate cuboids
  !>
  !> Requires
  !>   - mesh % n_elem
  !>   - mesh % element % face
  !>   - mesh % element % neighbor
  !>   - mesh % element % geometry % po
  !>   - mesh % element % geometry % x_e
  !>
  !> Generates
  !>   - mesh % element % geometry % x_c
  !>   - mesh % element % geometry % dx_m

  module subroutine BuildCuboids(mesh)
    class(Mesh_3D), intent(inout)  :: mesh !< mesh partition

    type(StandardElementOperators_1D), allocatable :: sop(:)
    real(RNP), allocatable, target :: work(:)
    real(RNP), pointer :: VL_inv(:,:)
    integer :: e, po

    allocate(sop( mesh%p_geom ), work( (mesh%p_geom + 1)**2 ))

    ! cuboids ..................................................................

    !$omp do
    do e = 1, mesh % n_elem
      associate( geometry => mesh % element(e) % geometry )

        po = geometry % po

        ! cuboid from truncated Legendre representation
        if (sop(po) % po /= po) then
          sop(po) = StandardElementOperators_1D(po)
        end if
        VL_inv(0:po,0:po) => work(1:(po+1)**2)
        call sop(po) % Get_Inverse_Legendre_VDM(VL_inv)
        call BuildElementCuboid(po, VL_inv, geometry % x_e, geometry % x_c)

      end associate
    end do

    ! mean spacing .............................................................

    call ComputeMeanSpacing(mesh)

  end subroutine BuildCuboids

  !-----------------------------------------------------------------------------
  !> Build cuboid from truncated Legendre series

  subroutine BuildElementCuboid(po, VL_inv, x, y)
    integer  , intent(in)  :: po                  !< polynomial order
    real(RNP), intent(in)  :: VL_inv(0:po,0:po)   !< inv Vandermonde matrix
    real(RNP), intent(in)  :: x(0:po,0:po,0:po,3) !< element Lobatto points
    real(RNP), intent(out) :: y(0:3,3)            !< cuboid coefficients

    real(RNP) :: a(0:po,0:po,3), b(0:po,3), c(0:1,0:1,0:1,3)
    integer   :: d, i, j, k, n

    n = po + 1

    do i = 0, 1
      a = reshape( matmul(VL_inv(i,:), reshape(x, [n,n*n*3])), [n,n,3] )
      do j = 0, 1
        b = reshape( matmul(VL_inv(j,:), reshape(a, [n,n*3])), [n,3] )
        do k = 0, 1
          c(i,j,k,1) = dot_product(VL_inv(k,:), b(:,1))
          c(i,j,k,2) = dot_product(VL_inv(k,:), b(:,2))
          c(i,j,k,3) = dot_product(VL_inv(k,:), b(:,3))
        end do
      end do
    end do

    do d = 1, 3
      y(0:3,d) = [ c(0,0,0,d), c(1,0,0,d), c(0,1,0,d), c(0,0,1,d) ]
    end do

  end subroutine BuildElementCuboid

  !-----------------------------------------------------------------------------
  !> Computes the mean spacing normal to faces

  subroutine ComputeMeanSpacing(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    type(ElementTransferBuffer_3D), allocatable, asynchronous, save :: buf_dx
    real(RNP), allocatable, save :: dx(:,:,:,:)
    integer :: e, f, i, n

    ! initialization ...........................................................

    !$omp master
    allocate(dx(3, 1, 1, mesh%n_elem + mesh%n_ghost))
    buf_dx = ElementTransferBuffer_3D(mesh, dx)
    !$omp end master
    !$omp barrier

    ! local cuboid spacing in ξ, η and ζ directions ............................

    !$omp do
    do e = 1, mesh % n_elem
      associate(x_c => mesh % element(e) % geometry % x_c)
        do i = 1, 3
          dx(i,1,1,e) = 2 * sqrt( x_c(i,1)**2 &
                                + x_c(i,2)**2 &
                                + x_c(i,3)**2 )
        end do
      end associate
    end do

    ! transfer cuboid spacing to ghosts ........................................

    call buf_dx % Transfer(mesh, dx, 1000)
    call buf_dx % Merge(dx)

    ! compute mean normal spacing ..............................................

    !$omp do
    do e = 1, mesh % n_elem
      associate( element => mesh % element(e) &
               , dx_m    => mesh % element(e) % geometry % dx_m )

        do f = 1, 6

          select case(f)
          case(1,2)
            dx_m(f) = dx(1,1,1,e)
          case(3,4)
            dx_m(f) = dx(2,1,1,e)
          case default
            dx_m(f) = dx(3,1,1,e)
          end select

          i = element % face(f) % i_neighbor

          if (i > 0) then

            n = element % neighbor(i) % id
            if (n == e) cycle

            select case(element % neighbor(i) % component) ! neighbor face
            case(1,2)
              dx_m(f) = 2 / (1/dx_m(f) + 1/dx(1,1,1,n))
            case(3,4)
              dx_m(f) = 2 / (1/dx_m(f) + 1/dx(2,1,1,n))
            case default
              dx_m(f) = 2 / (1/dx_m(f) + 1/dx(3,1,1,n))
            end select

          end if

        end do
      end associate
    end do

    ! clean-up .................................................................

    !$omp master
    deallocate(dx, buf_dx)
    !$omp end master

  end subroutine ComputeMeanSpacing

  !=============================================================================

end submodule MP_BuildCuboids

