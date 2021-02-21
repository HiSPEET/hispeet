!> summary:  Initialize elliptic operator for IP/DG-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2018/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build new elliptic operator for IP/DG-SEM
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Init
  use Array_Assignments, only: SetArray
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
!> Initialize operator with a combination of constant isotropic diffusivity and
!> constant isotropic spectral diffusivity

module subroutine Init_CI_svv(this, mesh, lambda, nu, nu_svv, bc, ip_opt,      &
  schwarz_opt)

! arguments ..................................................................

class(EllipticOperator3D_IP), intent(inout) :: this
class(MeshPartition), target, intent(in)    :: mesh   !< mesh partition
real(RNP),                    intent(in)    :: lambda !< Helmholtz parameter
real(RNP),                    intent(in)    :: nu     !< diffusivity
real(RNP),                    intent(in)    :: nu_svv !< SVV diffusivity
character,                    intent(in)    :: bc(:)  !< boundary conditions

!> options for the IP/DG method, including polynomial order `po` and `penalty`
class(IP_ElementOptions1D), intent(in) :: ip_opt

!> options for the Schwarz method
class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

! problem and discretization parameters ......................................

this % mesh      => mesh
this % lambda    =  lambda
this % nu_ci     =  nu
this % nu_ci_svv =  nu_svv
this % bc        =  bc

allocate(this % eop, source = IP_ElementOperators1D(ip_opt))

! Schwarz method .............................................................

if (present(schwarz_opt)) then
  allocate(this % schwarz)
  call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, nu_svv, bc)
end if

end subroutine Init_CI_svv

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

  type(TraceOperator), asynchronous, allocatable, save :: trace_op
  real(RNP), allocatable, save :: tr_nu(:,:,:,:)
  integer :: i, j, k, po

  ! problem and discretization parameters ......................................

  this % mesh   => mesh
  this % lambda =  lambda
  this % bc     =  bc

  allocate(this % eop, source = IP_ElementOperators1D(ip_opt))

  po = this % eop % po

  allocate(this % nu_vi(0:po,0:po,0:po,mesh%ne))
  call SetArray(this % nu_vi, nu)

  ! start generating traces of nu
  allocate(this % nu_hat(0:po,0:po,mesh%nf))
  allocate(tr_nu(0:po,0:po,2,mesh%nf))
  allocate(trace_op)
  call trace_op % GetTrace_Start(mesh, nu, tr_nu, tag=1000)

  ! Schwarz method .............................................................

  if (present(schwarz_opt)) then
    allocate(this % schwarz)
    call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, bc)
  end if

  ! max diffusivity on faces ...................................................

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

  ! clean-up ...................................................................

  deallocate(trace_op)
  deallocate(tr_nu)

end subroutine Init_VI

!===============================================================================

end submodule MP_Init
