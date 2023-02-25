!> summary:  Parent-to-child interpolation operators
!> author:   Joerg Stiller
!> date:     2023/01/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Parent_To_Child_Interpolation__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Standard_Operators__1D

  !-----------------------------------------------------------------------------
  !> Parent-to-child interpolation operator
  !>
  !> Provides the operators `A(0:po,0:po,c)` for interpolating parent data
  !> to children `c=1:2`.

  type ParentToChildInterpolation_1D
    real(RNP), allocatable :: A(:,:,:) !< interpolation operators
  end type ParentToChildInterpolation_1D

  ! constructor interface
  interface ParentToChildInterpolation_1D
    module procedure New_Points
    module procedure New_StdOps
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New parent-to-child interpolation from given basis type and points

  type(ParentToChildInterpolation_1D) function New_Points(basis, xi) result(this)
    character, intent(in)  :: basis  !< basis type
    real(RNP), intent(in)  :: xi(0:) !< parent basis points

    allocate(this % A(0:ubound(xi,1), 0:ubound(xi,1), 2))
    call BuildParentToChildInterpolation(basis, xi, this % A)

  end function New_Points

  !-----------------------------------------------------------------------------
  !> New parent-to-child interpolation from standard operators

  type(ParentToChildInterpolation_1D) function New_StdOps(eop) result(this)
    class(StandardOperators_1D), intent(in) :: eop !< standard operators

    allocate(this % A(0:eop%po, 0:eop%po, 2))
    call BuildParentToChildInterpolation(eop % basis, eop % x, this % A)

  end function New_StdOps

  !-----------------------------------------------------------------------------
  !> Build parent-to-child interpolation

  subroutine BuildParentToChildInterpolation(basis, xi, A)
    character, intent(in)  :: basis       !< basis type
    real(RNP), intent(in)  :: xi(0:)      !< parent basis points
    real(RNP), intent(out) :: A(0:,0:,1:) !< interpolation to child 1 and 2

    integer :: i, j, po

    po = ubound(xi, 1)

    select case(basis)
    case('G') ! Gauss
      do j = 0, po
      do i = 0, po
        A(i,j,1) = GaussPolynomial(j, xi, HALF * (xi(j) - ONE))
        A(i,j,2) = GaussPolynomial(j, xi, HALF * (xi(j) + ONE))
      end do
      end do
    case('R') ! Radau
      do j = 0, po
      do i = 0, po
        A(i,j,1) = RadauPolynomial(j, xi, HALF * (xi(j) - ONE))
        A(i,j,2) = RadauPolynomial(j, xi, HALF * (xi(j) + ONE))
      end do
      end do
    case default
      do j = 0, po
      do i = 0, po
        A(i,j,1) = LobattoPolynomial(j, xi, HALF * (xi(j) - ONE))
        A(i,j,2) = LobattoPolynomial(j, xi, HALF * (xi(j) + ONE))
      end do
      end do
    end select

  end subroutine BuildParentToChildInterpolation

  !=============================================================================

end module Parent_To_Child_Interpolation__1D
