!> summary:  Set problem parameters for elliptic operator with IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2019/02/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build new elliptic operator for IP/DG-SEM
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_SetProblem
  use CART__Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> Initialize problem with constant isotropic diffusivity

module subroutine SetProblem_CI(this, lambda, nu, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda !< Helmholtz parameter
  real(RNP), intent(in) :: nu     !< diffusivity
  character, intent(in) :: bc(:)  !< boundary conditions

  if (allocated(this % nu_vi )) deallocate(this % nu_vi )
  if (allocated(this % nu_hat)) deallocate(this % nu_hat)

  this % lambda    = lambda
  this % nu_ci     = nu
  this % bc        = bc

  if (allocated(this % schwarz)) then
    call this % schwarz % SetProblem(this%mesh, lambda, nu, bc)
  end if

end subroutine SetProblem_CI

!-------------------------------------------------------------------------------
!> Initialize problem with a combination of constant isotropic diffusivity and
!> constant isotropic spectral diffusivity

module subroutine SetProblem_CI_svv(this, lambda, nu, nu_svv, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda !< Helmholtz parameter
  real(RNP), intent(in) :: nu     !< diffusivity
  real(RNP), intent(in) :: nu_svv !< spectral diffusivity
  character, intent(in) :: bc(:)  !< boundary conditions

  if (allocated(this % nu_vi )) deallocate(this % nu_vi )
  if (allocated(this % nu_hat)) deallocate(this % nu_hat)

  this % lambda    = lambda
  this % nu_ci     = nu
  this % nu_ci_svv = nu_svv
  this % bc        = bc

  if (allocated(this % schwarz)) then
    write(*,*) "the Schwarz method does not work with SVV yet"
  end if

end subroutine SetProblem_CI_svv

!-------------------------------------------------------------------------------
!> Initialize problem with variable isotropic diffusivity

module subroutine SetProblem_VI(this, lambda, nu, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda         !< Helmholtz parameter
  real(RNP), intent(in) :: nu(0:,0:,0:,:) !< diffusivity
  character, intent(in) :: bc(:)          !< boundary conditions

  ! local variables ............................................................

  type(TraceOperator), allocatable, save :: trace_op
  real(RNP),           allocatable, save :: tr_nu(:,:,:,:)
  integer :: i, j, k, po, np

  associate(mesh => this % mesh)

    ! preparations ............................................................

    po = this % eop % po
    np = po + 1

    if (allocated(this % nu_ci))     deallocate(this % nu_ci)
    if (allocated(this % nu_ci_svv)) deallocate(this % nu_ci_svv)

    ! (re)allocate nu_vi, if necessary
    if (allocated(this % nu_vi)) then
      if (any(shape(this % nu_vi) /= shape(nu))) deallocate(this % nu_vi)
    end if
    if (.not. allocated(this % nu_vi)) then
      allocate(this % nu_vi, mold = nu)
    end if

    ! (re)allocate nu_hat, if necessary
    if (allocated(this % nu_hat)) then
      if (any(shape(this % nu_hat) /= [np,np,mesh%nf])) then
        deallocate(this % nu_hat)
      end if
    end if
    if (.not. allocated(this % nu_hat)) then
      allocate(this % nu_hat(0:po,0:po,mesh%nf))
    end if

    ! start generating traces of nu
    allocate(tr_nu(0:po,0:po,2,mesh%nf))
    allocate(trace_op)
    call trace_op % GetTrace_Start(mesh, nu, tr_nu, tag=1000)

    ! components ...............................................................

    this % lambda = lambda
    this % bc     = bc

    call SetArray(this % nu_vi, nu)

    if (allocated(this % schwarz)) then
      call this % schwarz % SetProblem(this%eop, mesh, lambda, nu, bc)
    end if

    call trace_op % GetTrace_Finish(mesh, tr_nu)

    associate(nu_hat => this % nu_hat)
      do k = 1, mesh%nf
        do j = 0, po
        do i = 0, po
          nu_hat(i,j,k) = max(tr_nu(i,j,1,k), tr_nu(i,j,2,k))
        end do
        end do
      end do
    end associate

   ! clean-up .................................................................

    deallocate(trace_op)
    deallocate(tr_nu)

  end associate

end subroutine SetProblem_VI

!===============================================================================

end submodule MP_SetProblem
