!> summary:   Conjugate gradients for elliptic equation with equidistant DG
!> author:    Joerg Stiller
!> date:      2017/07/06
!> license:t  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Conjugate gradients for elliptic equation with equidistant DG
!===============================================================================

module CART__DG_Elliptic_CI_Conj_Grad

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Array_Assignments, only: SetArray, MergeArrays, CalibrateArray
  use Array_Reductions,  only: ScalarProduct

  use XMPI

  use CART__Mesh_Partition
  use CART__DG_Element_Operators
  use CART__DG_Elliptic_CI_Operator
  use CART__DG_Elliptic_CI_Residual

  implicit none
  private

  public :: ConjugateGradients

contains

!-------------------------------------------------------------------------------
!> Conjugate gradient solver for elliptic equation with equidistant DG

subroutine ConjugateGradients( mesh, eop, lambda, nu, bc, u, f, i_max, r_red, &
                               r_max, d_min, ni )

  ! arguments ..................................................................

  class(MeshPartition),         intent(in) :: mesh  !< mesh partition
  class(DG_ElementOperators3D), intent(in) :: eop   !< DG element operators

  real(RNP), intent(in)    :: lambda           !< Helmholtz parameter
  real(RNP), intent(in)    :: nu               !< diffusivity
  character, intent(in)    :: bc(:)            !< boundary conditions {P,D,N}
  real(RNP), intent(inout) :: u(:,:,:,:)       !< approximate solution
  real(RNP), intent(in)    :: f(:,:,:,:)       !< right hand side
  integer,   intent(in)    :: i_max            !< maximum number of iterations
  real(RNP), intent(in)    :: r_red            !< targeted residual reduction
  real(RNP), optional, intent(in)  :: r_max    !< admissible residual
  real(RNP), optional, intent(in)  :: d_min    !< terminal threshold for Δr
  integer,   optional, intent(out) :: ni       !< number of iterations executed

  ! local variables ............................................................

  real(RNP), dimension(:,:,:,:), allocatable, save :: g, r, p, q
  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP) :: dd_min
  real(RNP) :: alpha, pq, rr, rr_old
  integer   :: i

  ! initialization .............................................................

  ! work space
  !$omp single
  allocate(g, mold = u)
  allocate(r, mold = u)
  allocate(p, mold = u)
  allocate(q, mold = u)
  !$omp end single

  !$acc data create(g,r,p,q) present(u,f)

  ! RHS
  call SetArray(g, f)

  ! calibrate RHS of singular problem
  if (abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')) then
    call CalibrateArray(g, mesh%comm)
  end if

  ! terminal threshold
  if (present(d_min)) then
    dd_min = d_min**2
  else
    dd_min = ZERO
  end if

  ! initial residual ...........................................................

  call EllipticResidual(mesh, eop, lambda, nu, bc, u, g, r)
  call SetArray(p, r)

  rr = ScalarProduct(r, r, mesh%comm)

  !$omp single
  rr_term  = max(ZERO, sqrt(rr) * r_red)**2
  if (present(r_max)) then
    rr_term = max(rr_term, max(ZERO, r_max)**2)
  end if
  !$omp end single

  rr_old = 0

  ! iteration ..................................................................

  do i = 1, i_max

    ! termination check  . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    ! MPI master decides about termination
    !$omp master
    if (mesh%part == 0) then
      converged = rr <= rr_term .or. abs(rr - rr_old) <= dd_min
    end if
    rr_old = rr
    call XMPI_Bcast(converged, root=0, comm=mesh%comm)
    !$omp end master
    !$omp barrier

    if (converged) exit

    ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    call EllipticOperator(mesh, eop, lambda, nu, bc, p, q)

    pq = ScalarProduct(p, q, mesh%comm)
    alpha = rr_old / pq

    call MergeArrays(ONE, u,  alpha, p)

    if (mod(i,50) == 0) then
      ! compute true residual to get rid of round-off errors
      call EllipticResidual(mesh, eop, lambda, nu, bc, u, g, r)
    else
      call MergeArrays(ONE, r, -alpha, q)
    end if

    rr = ScalarProduct(r, r, mesh%comm)

    call  MergeArrays(rr/rr_old, p, ONE, r)

  end do

  if (present(ni)) ni = i - 1

  !$acc end data

  !$omp barrier
  !$omp master
  deallocate(g, r, p, q)
  !$omp end master


end subroutine ConjugateGradients

!===============================================================================

end module CART__DG_Elliptic_CI_Conj_Grad
