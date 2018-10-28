!> summary:  Utilities for continuous cuboidal elements
!> author:   Joerg Stiller
!> date:     2018/10/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Utilities for continuous cuboidal elements
!>
!> Provides common routines for GLL-based nodal CG-SEM, including
!>
!>   *  computation of the global mass matrix (`GetMassMatrix`)  -- not yet
!>   *  averaging of discontinuous data (`MakeContinuous`)       -- not yet
!>   *  provision of point weights based on valency (`GetPointWeights`)
!>
!> Some routines require the specification of boundary condition types, which
!> are passed in the character array `bc(1:2)`. The following types are
!> supported:
!>
!>   *  periodic  (`'P'`)
!>   *  Dirichlet (`'D'`)
!>   *  Neumann   (`'N'`)
!>
!===============================================================================

module CART__CG_Element_Utilities
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error
  use Array_Assignments
  use CART__Mesh_Partition
  use CART__Element_Transfer_Buffer
  implicit none
  private

  public :: GetPointWeights

contains

!-------------------------------------------------------------------------------
!> Generates point weights based on valency and boundary conditions
!>
!> Computes the weights required for evaluating dot products in continuous
!> spectral element methods. In general, the weights are just the inverse of
!> the point valency. If the BC types are given, the weights are set to zero
!> for all Dirichlet points.

subroutine GetPointWeights(mesh, bc, w)
  class(MeshPartition), intent(in)  :: mesh          !< mesh parition
  character,  optional, intent(in)  :: bc(:)         !< BC types
  real(RNP),            intent(out) :: w(0:,0:,0:,:) !< point weights

  integer, allocatable :: v(:,:,:,:)
  integer :: b, e, f, i, j, k, l, po

  ! initialization .............................................................

  po = ubound(v, 1)

  allocate(v(0:po,0:po,0:po,mesh%ne))
  call mesh % GetPointValency('GLL', v)

  ! weights based on valency ...................................................

  do e = 1, mesh%ne
    do k = 0, po
    do j = 0, po
    do i = 0, po
      w(i,j,k,e) = ONE / v(i,j,k,e)
    end do
    end do
    end do
  end do

  ! zero Dirichlet weights .....................................................

  if (present(bc)) then

    Boundaries: do b = 1, size(bc)

      Boundary_Faces: associate(face => mesh%boundary(b)%face)

        if (bc(b) /= 'D') cycle

        do l = 1, size(face)
          e = face(l) % mesh_element % id
          f = face(l) % mesh_element % face

          select case(f)

          case(1) ! direction 1: west
            do k = 0, po
            do j = 0, po
              w(0,j,k,e) = 0
            end do
            end do

          case(2) ! direction 1: east
            do k = 0, po
            do j = 0, po
              w(po,j,k,e) = 0
            end do
            end do

          case(3) ! direction 2: south
            do k = 0, po
            do i = 0, po
              w(i,0,k,e) = 0
            end do
            end do

          case(4) ! direction 2: north
            do k = 0, po
            do i = 0, po
              w(i,po,k,e) = 0
            end do
            end do

          case(5) ! direction 3: bottom
            do j = 0, po
            do i = 0, po
              w(i,j,0,e) = 0
            end do
            end do

          case(6) ! direction 3: top
            do j = 0, po
            do i = 0, po
              w(i,j,po,e) = 0
            end do
            end do

          end select
        end do

      end associate Boundary_Faces
    end do Boundaries

  end if

end subroutine GetPointWeights

!===============================================================================

end module CART__CG_Element_Utilities
