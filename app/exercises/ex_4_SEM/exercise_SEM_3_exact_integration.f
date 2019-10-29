!> \file       exercise_3_exact_integration.f 
!> \brief      Test program examining the integration error on GLL points
!> \author     Immo Huismann
!> \date       2013/12/02
!> \copyright  Institute of Fluid Mechanics &
!>             TU Dresden, 01062 Dresden, Germany
!> \details
!> This program integrates the product of two polynomials of degree p and
!> calculates the solution error. The approximative and the exact mass matrix
!> are to be used for the integration, but are not implemented correctly.
!> For this program to work the TODOs 1 and 2 have to be done.
!===============================================================================

program Exercise_3_Exact_Integration
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF, TWO
  use Gauss_Jacobi, &
       only: GLL_Points, GLL_Weights, GL_Points, GL_Weights, GLL_Polynomial
  implicit none

  integer,   parameter   :: PO_MAX = 24 ! maximum polynomial degree to test

  real(RNP), allocatable :: xi(:)                  ! GLL points
  real(RNP), allocatable :: u_h(:)                 ! interpolation coefficients

  real(RNP), allocatable :: M_ex(:,:), M_appr(:,:) ! mass matrices

  real(RNP) :: int_appr(PO_MAX), int_ex(PO_MAX)    ! integrals
  real(RNP) :: err_appr(PO_MAX), err_ex(PO_MAX)    ! errors

  integer :: po   ! the polynomial degree for the collocation points
  integer :: uni  ! output unit

! Error calculation ............................................................

  ! check the integration error on all polynomial degrees < PO_MAX
  do po = 1, PO_MAX

    allocate(xi(0:po), M_ex(0:po,0:po), M_appr(0:po,0:po), u_h(0:po))

    ! get the GLL points
    xi = GLL_Points (po)

    ! set the mass matrices
    M_appr = MassMatrixApproximative(po)
    M_ex   = MassMatrixExact        (po)

    ! interpolate a monomial of degree po
    u_h = Monomial(po, xi)

    ! integrate u_h * u_h
    int_appr(po) = dot_product(u_h,matmul(M_appr,u_h))
    int_ex  (po) = dot_product(u_h,matmul(M_ex  ,u_h))
    err_appr(po) = abs(int_appr(po) - 1)
    err_ex  (po) = abs(int_ex  (po) - 1)

    deallocate(xi, M_ex, M_appr, u_h)

  end do

  ! Output .....................................................................
  ! print the error into a file
  open(newunit=uni, file='output_exercise_3_exact_integration.dat')
  do po = 1, PO_MAX
    write(uni,'(I4,4E15.7)') po,int_appr(po),int_ex(po),err_appr(po),err_ex(po)
  end do
  close(uni)

contains

!-------------------------------------------------------------------------------
!> \brief   A monomial whose squares integral on the standard element is one.
!> \author  Immo Huismann

function Monomial(p, x) result(y)
  integer  , intent(in) :: p              !> Degree of the monomial
  real(RNP), intent(in) :: x(:)           !> Node locations
  real(RNP)             :: y(ubound(x,1)) !> Values on node locations

  y = sqrt(HALF * (2 * p + 1)) * x ** p

end function Monomial

!-------------------------------------------------------------------------------
!> \brief   Calculates the approximative mass matrix for the GLL points
!> \author  Immo Huismann

function MassMatrixApproximative(po) result(M)
  integer,  intent(in) :: po           !> The polynomial degree
  real(RNP)            :: M(0:po,0:po) !> The mass matrix to set

  real(RNP)            :: xi(0:po) ! GLL nodes
  real(RNP)            :: wi(0:po) ! GLL Weights
  integer              :: i        ! Loop index
  
  xi = GLL_Points (po)
  wi = GLL_Weights(xi)

  ! TODO 1 - Set the approximative mass matrix
  
  ! Current implementation: Weights from a first-order spline
  M = 0
  do i = 0, po
    M(i,i) = TWO / (po + 1)
  end do

end function MassMatrixApproximative

!-------------------------------------------------------------------------------
!> \brief   Calculates the exact mass matrix for the GLL points
!> \author  Immo Huismann

pure function MassMatrixExact(po) result(M)
  integer,  intent(in) :: po           !> Polynomial degree
  real(RNP)            :: M(0:po,0:po) !> Mass matrix to set

  real(RNP)            :: xi_GLL(0:po) ! GLL Nodes
  real(RNP)            :: xi_GL (0:po) ! GL Nodes
  real(RNP)            :: wi_GL (0:po) ! GL Weights
  integer              :: i, j         ! loop indices
  real(RNP)            :: S(0:po,0:po) ! interpolation matrix
  
  ! set the GLL and GL variables
  xi_GLL = GLL_Points(po)
  xi_GL  = GL_Points (po)
  wi_GL  = GL_Weights(xi_GL)

  ! TODO 2 - set the interpolation matrix and the mass matrix correctly

  ! set the interpolation matrix
  S = 0

  do j = 0, po
  do i = 0, po
    if(i == j) S(i,i) = 1
  end do
  end do

  ! set M
  M = 0
  do i = 0, po
    M(i,i) = wi_GL(i)
  end do

end function MassMatrixExact

end program Exercise_3_Exact_Integration
