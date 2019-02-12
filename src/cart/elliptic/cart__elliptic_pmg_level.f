!> summary:  Polynomial multigrid level for use with elliptic solvers
!> author:   Joerg Stiller
!> date:     2019/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid level for use with DG elliptic solvers
!===============================================================================

module CART__Elliptic_PMG_Level
  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  use TPO_AAA
  use Standard_Operators_1D
  use IP_Element_Operators_1D
  use CART__Mesh_Partition
  use CART__Schwarz_Operator
  use CART__Elliptic_Operator
  use CART__Elliptic_Operator_IP

  implicit none
  private

  public :: PMG_Level

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid level

  type PMG_Level

    integer :: po = -1  !< polynomial order
    integer :: ne = -1  !< number of local elements
    integer :: ns1 = 1  !< number of pre-smoothing steps
    integer :: ns2 = 1  !< number of post-smoothing steps

    class(EllipticOperator3D), allocatable :: elliptic_op !< elliptic operator

    real(RNP), allocatable :: u(:,:,:,:) !< solution
    real(RNP), allocatable :: f(:,:,:,:) !< RHS
    real(RNP), allocatable :: v(:,:,:,:) !< residual / correction

    real(RNP), allocatable, private :: p2f_op(:,:) ! prolongation operator
    real(RNP), allocatable, private :: i2c_op(:,:) ! interpolation operator
    real(RNP), allocatable, private :: r2c_op(:,:) ! restriction operator
    real(RNP), allocatable, private :: t2c_op(:,:) ! truncation operator

  contains

    generic :: Init_TopLevel => Init_TopLevel__IP_CI, Init_TopLevel__IP_VI
    procedure, private :: Init_TopLevel__IP_CI
    procedure, private :: Init_TopLevel__IP_VI

    procedure :: Init_CoarseLevel

    procedure :: Prolongate  => C2F_Prolongation
    procedure :: Interpolate => F2C_Interpolation
    procedure :: Restrict    => F2C_Restriction
    procedure :: Truncate    => F2C_Truncation

  end type PMG_Level

contains

!===============================================================================
! Initialization of PMG_Level

!-------------------------------------------------------------------------------
!> Initialization of top level with IP/DG-SEM and constant isotropic diffusivity

subroutine Init_TopLevel__IP_CI( this, ns1, ns2, mesh, lambda, nu, bc, &
                                 ip_opt, schwarz_opt, po_parent        )

  ! arguments ..................................................................

  class(PMG_Level),     intent(inout) :: this
  integer,              intent(in)    :: ns1    !< num pre-smoothing steps
  integer,              intent(in)    :: ns2    !< num post-smoothing steps
  class(MeshPartition), intent(in)    :: mesh   !< mesh partition
  real(RNP),            intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),            intent(in)    :: nu     !< diffusivity
  character,            intent(in)    :: bc(:)  !< boundary conditions

  class(IP_ElementOptions1D), intent(in) :: ip_opt
  !< options for the IP/DG method, including polynomial order `po` and `penalty`

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  integer, optional, intent(in) :: po_parent
  !< order of next coarser (parent) level, if any

  ! parameters .................................................................

  this % po  = ip_opt % po
  this % ne  = mesh % ne
  this % ns1 = ns1
  this % ns2 = ns2

  ! elliptic operators .........................................................

  this % elliptic_op = EllipticOperator3D_IP( mesh, lambda, nu, bc, &
                                              ip_opt, schwarz_opt   )

  ! transfer operators .........................................................

  if (present(po_parent)) then
    call Build_F2C_TransferOps(this, po_parent)
  end if

  ! workspace ..................................................................

  call GetWorkspace(this)

end subroutine Init_TopLevel__IP_CI

!-------------------------------------------------------------------------------
!> Initialization of top level with IP/DG-SEM and variable isotropic diffusivity

subroutine Init_TopLevel__IP_VI( this, ns1, ns2, mesh, lambda, nu, bc, &
                                 ip_opt, schwarz_opt, po_parent        )

  ! arguments ..................................................................

  class(PMG_Level),     intent(inout) :: this
  integer,              intent(in)    :: ns1            !< n pre-smoothing steps
  integer,              intent(in)    :: ns2            !< n post-smoothing steps
  class(MeshPartition), intent(in)    :: mesh           !< mesh partition
  real(RNP),            intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP),            intent(in)    :: nu(0:,0:,0:,:) !< diffusivity
  character,            intent(in)    :: bc(:)          !< boundary conditions

  class(IP_ElementOptions1D), intent(in) :: ip_opt
  !< options for the IP/DG method, including polynomial order `po` and `penalty`

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  integer, optional, intent(in) :: po_parent
  !< order of next coarser (parent) level, if any

  ! parameters .................................................................

  this % po  = ip_opt % po
  this % ne  = mesh % ne
  this % ns1 = ns1
  this % ns2 = ns2

  ! elliptic operators .........................................................

  this % elliptic_op = EllipticOperator3D_IP( mesh, lambda, nu, bc, &
                                              ip_opt, schwarz_opt   )

  ! transfer operators .........................................................

  if (present(po_parent)) then
    call Build_F2C_TransferOps(this, po_parent)
  end if

  ! workspace ..................................................................

  call GetWorkspace(this)

end subroutine Init_TopLevel__IP_VI

!-------------------------------------------------------------------------------
!> Initialization of coarse level

subroutine Init_CoarseLevel(this, po, ns1, ns2, child, schwarz_opt, po_parent)

  ! arguments ..................................................................

  class(PMG_Level), intent(inout) :: this
  integer,          intent(in)    :: po    !< polynomial order
  integer,          intent(in)    :: ns1   !< num pre-smoothing steps
  integer,          intent(in)    :: ns2   !< num post-smoothing steps
  class(PMG_Level), intent(in)    :: child !< next finer (child) level

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  integer, optional, intent(in) :: po_parent
  !< order of next coarser (parent) level, if any

  ! internal variables .........................................................

  real(RNP), allocatable, save :: nu_vi(:,:,:,:)

  ! parameters .................................................................

  this % po  = po
  this % ne  = child % ne
  this % ns1 = ns1
  this % ns2 = ns2

  ! elliptic operators .........................................................

  if (allocated(child % elliptic_op % nu_vi)) then
    allocate(nu_vi(0:this%po, 0:this%po, 0:this%po, this%ne))
  end if

  select type(elliptic_op => child % elliptic_op)

  class is(EllipticOperator3D_IP)

    block
      type(IP_ElementOptions1D) :: ip_opt

      ip_opt = IP_ElementOptions1D(elliptic_op % eop, po)

      if (allocated(child % elliptic_op % nu_ci)) then

        this % elliptic_op = EllipticOperator3D_IP( elliptic_op % mesh,   &
                                                    elliptic_op % lambda, &
                                                    elliptic_op % nu_ci,  &
                                                    elliptic_op % bc,     &
                                                    ip_opt,               &
                                                    schwarz_opt           )

      else if (allocated(nu_vi)) then

        !call child % Truncate(elliptic_op % nu_vi, nu_vi)
        call child % Interpolate(elliptic_op % nu_vi, nu_vi)

        this % elliptic_op = EllipticOperator3D_IP( elliptic_op % mesh,   &
                                                    elliptic_op % lambda, &
                                                    nu_vi,                &
                                                    elliptic_op % bc,     &
                                                    ip_opt,               &
                                                    schwarz_opt           )
      end if
    end block
  end select

  if (allocated(nu_vi)) then
    deallocate(nu_vi)
  end if

  ! transfer operators .........................................................

  call Build_C2F_TransferOps(this, child)

  if (present(po_parent)) then
    call Build_F2C_TransferOps(this, po_parent)
  end if

  ! workspace ..................................................................

  call GetWorkspace(this)

end subroutine Init_CoarseLevel

!-------------------------------------------------------------------------------
!> Prolongation: coarse-to-fine interpolation

subroutine C2F_Prolongation(this, uc, uf)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uc(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uf(0:,0:,0:,:) !< fine (child) mesh variable

  call TPO_AAA_Eval(size(uf,1), size(uc,1), size(uc,4), this%p2f_op, uc, uf)

end subroutine C2F_Prolongation

!-------------------------------------------------------------------------------
!> Fine-to-coarse interpolation

subroutine F2C_Interpolation(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA_Eval(size(uc,1), size(uf,1), size(uf,4), this%i2c_op, uf, uc)

end subroutine F2C_Interpolation

!-------------------------------------------------------------------------------
!> Fine-to-coarse restriction -- no weighting / no assembly !!

subroutine F2C_Restriction(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA_Eval(size(uc,1), size(uf,1), size(uf,4), this%r2c_op, uf, uc)

end subroutine F2C_Restriction

!-------------------------------------------------------------------------------
!> Fine-to-coarse restriction -- no weighting / no assembly !!

subroutine F2C_Truncation(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA_Eval(size(uc,1), size(uf,1), size(uf,4), this%t2c_op, uf, uc)

end subroutine F2C_Truncation

!===============================================================================
! Helpers

!-------------------------------------------------------------------------------
!> Build coarse-to-fine transfer operators

subroutine Build_C2F_TransferOps(this, fine)
  class(PMG_Level), intent(inout) :: this !< current (coarse) level
  class(PMG_Level), intent(in)    :: fine !< fine level

  real(RNP) :: xc(0:this%po), xf(0:fine%po)
  integer   :: i, k, pc, pf

  pc = this % po
  pf = fine % po
  xc = this % elliptic_op % eop % x
  xf = fine % elliptic_op % eop % x

  ! prolongation to fine level ...............................................

  allocate(this % p2f_op(0:pf, 0:pc))
  do k = 0, pc
  do i = 0, pf
    this % p2f_op(i,k) = GLL_Polynomial(k, xc, xf(i))
  end do
  end do

end subroutine Build_C2F_TransferOps

!-------------------------------------------------------------------------------
!> Build fine-to-coarse transfer operators

subroutine Build_F2C_TransferOps(this, pc)
  class(PMG_Level), intent(inout) :: this !< current (fine) level
  integer,          intent(in)    :: pc   !< polynomial order of coarse level

  real(RNP), allocatable :: VL(:,:), VL_inv(:,:)
  real(RNP) :: xf(0:this%po), xc(0:pc)
  integer   :: i, k, pf

  pf = this % po
  xf = this % elliptic_op % eop % x
  xc = GLL_Points(pc)

  ! interpolation to child level ...............................................

  allocate(this % i2c_op(0:pc,0:pf))
  do k = 0, pf
  do i = 0, pc
    this % i2c_op(i,k) = GLL_Polynomial(k, xf, xc(i))
  end do
  end do

  ! restriction to child level .................................................

  allocate(this % r2c_op(0:pc,0:pf))
  do k = 0, pf
  do i = 0, pc
    this % r2c_op(i,k) = GLL_Polynomial(i, xc, xf(k))
  end do
  end do

  ! truncation to child level ..................................................

  ! initialize with fine-to-coarse interpolation
  allocate(this % t2c_op(0:pc,0:pf), source = this % i2c_op)

  ! prepend order reduction, if possible
  associate(eop => this % elliptic_op % eop, t2c_op => this % t2c_op)
    if (eop % HasLegendreVDM()) then
      allocate(VL(0:pf,0:pf), VL_inv(0:pf,0:pf))
      call eop % GetLegendreVDM(VL)
      call eop % GetInverseLegendreVDM(VL_inv)
      t2c_op = matmul(t2c_op, matmul(VL(:,0:pc), VL_inv(0:pc,:)))
    end if
  end associate

end subroutine Build_F2C_TransferOps

!-------------------------------------------------------------------------------
!> Allocate workspace

subroutine GetWorkspace(this)
  class(PMG_Level), intent(inout) :: this

  associate(po => this%po, ne => this%ne)

    if (allocated(this % u)) then
      if (any(shape(this % u) /= [po+1, po+1, po+1, ne])) then
        deallocate(this % u)
        deallocate(this % f)
        deallocate(this % v)
      end if
    end if

    if (.not. allocated(this % u)) then
      allocate(this % u(0:po, 0:po, 0:po, ne))
      allocate(this % f, mold = this % u)
      allocate(this % v, mold = this % u)
    end if

  end associate

end subroutine GetWorkspace

!===============================================================================

end module CART__Elliptic_PMG_Level
