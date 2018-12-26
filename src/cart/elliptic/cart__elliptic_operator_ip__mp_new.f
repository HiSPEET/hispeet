!> summary:  Build new elliptic operator for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2018/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build new elliptic operator for IP/DG-SEM
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_New
  use Array_Assignments, only: AssignArray
  use CART__Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> New operator with constant isotropic diffusivity

module subroutine New_CI(this, mesh, lambda, nu, bc, po, penalty, schwarz_opt)

  ! arguments ..................................................................

  class(EllipticOperator3D_IP), intent(inout) :: this

  class(MeshPartition), target, intent(in) :: mesh        !< mesh partition
  real(RNP),                    intent(in) :: lambda      !< Helmholtz parameter
  real(RNP),                    intent(in) :: nu          !< diffusivity
  character,                    intent(in) :: bc(:)       !< boundary conditions
  integer,                      intent(in) :: po          !< polynomial order
  real(RNP),                    intent(in) :: penalty     !< penalty parameter
  class(SchwarzOptions3D),      intent(in) :: schwarz_opt !< Schwarz options

  optional :: penalty, schwarz_opt

  ! problem and discretization parameters ......................................

  this % mesh => mesh

  this % lambda = lambda
  this % nu_ci  = nu
  this % bc     = bc

  call this % eop % New(po, penalty)

  ! Schwarz method .............................................................

  if (present(schwarz_opt)) then
    allocate(this % schwarz)
    call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, bc)
  end if

end subroutine New_CI

!-------------------------------------------------------------------------------
!> New operator with variable isotropic diffusivity

module subroutine New_VI(this, mesh, lambda, nu, bc, po, penalty, schwarz_opt)

  ! arguments ..................................................................

  class(EllipticOperator3D_IP), intent(inout) :: this

  class(MeshPartition),    target, &
                           intent(in) :: mesh           !< mesh partition
  real(RNP),               intent(in) :: lambda         !< Helmholtz parameter
  real(RNP),               intent(in) :: nu(0:,0:,0:,:) !< diffusivity
  character,               intent(in) :: bc(:)          !< boundary conditions
  integer,                 intent(in) :: po             !< polynomial order
  real(RNP),               intent(in) :: penalty        !< penalty parameter
  class(SchwarzOptions3D), intent(in) :: schwarz_opt    !< Schwarz options

  optional :: penalty, schwarz_opt

  ! local variables ............................................................

  type(TraceOperator), allocatable, save :: trace_op
  real(RNP),           allocatable, save :: tr_nu(:,:,:,:)
  integer :: i, j, k

  ! problem and discretization parameters ......................................

  this % mesh => mesh

  this % lambda = lambda
  this % bc     = bc

  allocate(this % nu_vi, mold = nu)
  call AssignArray(this % nu_vi, nu)

  ! start generating traces of nu
  allocate(this % nu_f(0:po,0:po,mesh%nf))
  allocate(tr_nu(0:po,0:po,2,mesh%nf))
  allocate(trace_op)
  call trace_op % GetTrace_Start(mesh, nu, tr_nu, tag=1000)

  ! 1D IP/DG operators
  call this % eop % New(po, penalty)

  ! Schwarz method .............................................................

  if (present(schwarz_opt)) then
    allocate(this % schwarz)
    call this % schwarz % New(schwarz_opt, this%eop, mesh, lambda, nu, bc)
  end if

  ! max diffusivity on faces ...................................................

  call trace_op % GetTrace_Finish(mesh, tr_nu)

  associate(nu_f => this % nu_f)
    do k = 1, mesh%nf
      do j = 0, po
      do i = 0, po
        nu_f(i,j,k) = max(tr_nu(i,j,1,k), tr_nu(i,j,2,k))
      end do
      end do
    end do
  end associate

  ! clean-up ...................................................................

  deallocate(trace_op)
  deallocate(tr_nu)

end subroutine New_VI

!===============================================================================

end submodule MP_New
