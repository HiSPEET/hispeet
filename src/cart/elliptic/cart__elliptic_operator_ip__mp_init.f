!> summary:  Initialize elliptic operator for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2018/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build new elliptic operator for IP/DG-SEM
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Init
  use Array_Assignments, only: AssignArray
  use CART__Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> Initialize operator with constant isotropic diffusivity

module subroutine Init_CI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)

  ! arguments ..................................................................

  class(EllipticOperator3D_IP), intent(inout) :: this
  class(MeshPartition), target, intent(in)    :: mesh   !< mesh partition
  real(RNP),                    intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),                    intent(in)    :: nu     !< diffusivity
  character,                    intent(in)    :: bc(:)  !< boundary conditions

  !> options for the IP/DG method, including polynomial order `po` and `penalty`
  class(IP_ElementOptions1D), intent(in) :: ip_opt

  !> options for the Schwarz method
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  ! problem and discretization parameters ......................................

  this % mesh   => mesh
  this % lambda =  lambda
  this % nu_ci  =  nu
  this % bc     =  bc

  allocate(this % eop, source = IP_ElementOperators1D(ip_opt))

  ! Schwarz method .............................................................

  if (present(schwarz_opt)) then
    allocate(this % schwarz)
    call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, bc)
  end if

end subroutine Init_CI

!-------------------------------------------------------------------------------
!> Initialize operator with variable isotropic diffusivity

module subroutine Init_VI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)

  ! arguments ..................................................................

  class(EllipticOperator3D_IP), intent(inout) :: this
  class(MeshPartition), target, intent(in)    :: mesh  !< mesh partition
  real(RNP), intent(in) :: lambda                      !< Helmholtz parameter
  real(RNP), intent(in) :: nu(0:,0:,0:,:)              !< diffusivity
  character, intent(in) :: bc(:)                       !< boundary conditions

  !> options for the IP/DG method, including polynomial order `po` and `penalty`
  class(IP_ElementOptions1D), intent(in) :: ip_opt

  !> options for the Schwarz method
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  ! local variables ............................................................

  type(TraceOperator), allocatable, save :: trace_op
  real(RNP),           allocatable, save :: tr_nu(:,:,:,:)
  integer :: i, j, k, po

  ! problem and discretization parameters ......................................
!### CHECK
print *, '? 00'
!### CHECK END

  this % mesh   => mesh
  this % lambda =  lambda
  this % bc     =  bc

  allocate(this % eop, source = IP_ElementOperators1D(ip_opt))
!### CHECK
print *, '? 01'
!### CHECK END

  po = this % eop % po

  allocate(this % nu_vi, mold = nu)
  call AssignArray(this % nu_vi, nu)
!### CHECK
print *, '? 02'
!### CHECK END

  ! start generating traces of nu
  allocate(this % nu_hat(0:po,0:po,mesh%nf))
  allocate(tr_nu(0:po,0:po,2,mesh%nf))
  allocate(trace_op)
!### CHECK
print *, '? shape(nu)     =', shape(nu)
print *, '? shape(nu_hat) =', shape(this % nu_hat)
print *, '? shape(tr_nu)  =', shape(tr_nu)
!### CHECK END
  call trace_op % GetTrace_Start(mesh, nu, tr_nu, tag=1000)
!### CHECK
print *, '? 03'
!### CHECK END

  ! Schwarz method .............................................................

  if (present(schwarz_opt)) then
    allocate(this % schwarz)
    call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, bc)
  end if
!### CHECK
print *, '? 04'
!### CHECK END

  ! max diffusivity on faces ...................................................

  call trace_op % GetTrace_Finish(mesh, tr_nu)
!### CHECK
print *, '? 05'
!### CHECK END

  associate(nu_hat => this % nu_hat)
    do k = 1, mesh%nf
      do j = 0, po
      do i = 0, po
        nu_hat(i,j,k) = max(tr_nu(i,j,1,k), tr_nu(i,j,2,k))
      end do
      end do
    end do
  end associate
!### CHECK
print *, '? 06'
!### CHECK END

  ! clean-up ...................................................................

  deallocate(trace_op)
  deallocate(tr_nu)
!### CHECK
print *, '? 0X'
!### CHECK END

end subroutine Init_VI

!===============================================================================

end submodule MP_Init
