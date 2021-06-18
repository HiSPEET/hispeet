!> summary:  Polynomial multigrid level for use with elliptic solvers
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2019/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__Elliptic_PMG_Level
  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  use TPO__AAA__3D
  use Standard_Operators__1D
  use IP_Element_Operators__1D
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

    integer   :: po = -1        !< polynomial order
    integer   :: ne = -1        !< number of local elements
    integer   :: ns1 = 1        !< number of pre-smoothing steps
    integer   :: ns2 = 1        !< number of post-smoothing steps
    character :: smoother = 'S' !< 'S': Schwarz, 'C': CG, 'P': Schwarz-PCG

    class(EllipticOperator3D), allocatable :: elliptic_op !< elliptic operator

    real(RNP), allocatable :: u(:,:,:,:) !< solution
    real(RNP), allocatable :: f(:,:,:,:) !< RHS
    real(RNP), allocatable :: v(:,:,:,:) !< residual / correction

    real(RNP), allocatable, private :: p2f_op(:,:) ! prolongation operator
    real(RNP), allocatable, private :: i2c_op(:,:) ! interpolation operator
    real(RNP), allocatable, private :: r2c_op(:,:) ! restriction operator
    real(RNP), allocatable, private :: t2c_op(:,:) ! truncation operator

  contains

    procedure :: Init_TopLevel_IP
    procedure :: Init_CoarseLevel

    procedure :: SetTopLevelProblem_CI
    procedure :: SetTopLevelProblem_CI_svv
    procedure :: SetTopLevelProblem_VI
    procedure :: SetCoarseProblem

    procedure :: GetWorkspace
    procedure :: FreeWorkspace

    procedure :: Prolongate
    procedure :: Interpolate
    procedure :: Restrict
    procedure :: Truncate

  end type PMG_Level

contains

!===============================================================================
! Initialization of PMG_Level

!-------------------------------------------------------------------------------
!> Initialization of top level with IP/DG-SEM and no problem data

subroutine Init_TopLevel_IP( this, ns1, ns2, smoother, mesh, &
                             ip_opt , schwarz_opt, po_coarse )

  ! arguments ..................................................................

  class(PMG_Level),     intent(inout) :: this
  integer,              intent(in)    :: ns1      !< num pre-smoothing steps
  integer,              intent(in)    :: ns2      !< num post-smoothing steps
  character,            intent(in)    :: smoother !< smoothing method
  class(MeshPartition), intent(in)    :: mesh     !< mesh partition

  class(IP_ElementOptions_1D), intent(in) :: ip_opt
  !< options for the IP/DG method, including polynomial order `po` and `penalty`

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  integer, optional, intent(in) :: po_coarse
  !< order of next coarser level, if any

  ! parameters .................................................................

  this % po       = ip_opt % po
  this % ne       = mesh % ne
  this % ns1      = ns1
  this % ns2      = ns2
  this % smoother = smoother

  ! elliptic operators .........................................................

  this % elliptic_op = EllipticOperator3D_IP(mesh, ip_opt, schwarz_opt)

  ! transfer operators .........................................................

  if (present(po_coarse)) then
    call Build_F2C_TransferOps(this, po_coarse)
  end if

end subroutine Init_TopLevel_IP

!-------------------------------------------------------------------------------
!> Initialization of coarse level

subroutine Init_CoarseLevel(this, po, ns1, ns2, fine, schwarz_opt, po_coarse)

  ! arguments ..................................................................

  class(PMG_Level), intent(inout) :: this
  integer,          intent(in)    :: po    !< polynomial order
  integer,          intent(in)    :: ns1   !< num pre-smoothing steps
  integer,          intent(in)    :: ns2   !< num post-smoothing steps
  class(PMG_Level), intent(in)    :: fine  !< next finer level

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  integer, optional, intent(in) :: po_coarse
  !< order of next coarser (parent) level, if any

  ! parameters .................................................................

  this % po       = po
  this % ne       = fine % ne
  this % ns1      = ns1
  this % ns2      = ns2
  this % smoother = fine % smoother

  ! elliptic operators .........................................................

  select type(elliptic_op => fine % elliptic_op)

  class is(EllipticOperator3D_IP)

     this % elliptic_op = EllipticOperator3D_IP(                        &
                              elliptic_op % mesh,                       &
                              IP_ElementOptions_1D(elliptic_op%eop, po), &
                              schwarz_opt                               )
  end select

  ! transfer operators .........................................................

  call Build_C2F_TransferOps(this, fine)

  if (present(po_coarse)) then
    call Build_F2C_TransferOps(this, po_coarse)
  end if

end subroutine Init_CoarseLevel

!===============================================================================
! (Re)initialization of problem data

!-------------------------------------------------------------------------------
!> Initialize top-level problem with constant isotropic diffusivity

subroutine SetTopLevelProblem_CI(this, lambda, nu, bc)
  class(PMG_Level), intent(inout) :: this
  real(RNP),        intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),        intent(in)    :: nu     !< diffusivity
  character,        intent(in)    :: bc(:)  !< boundary conditions

  call this % elliptic_op % SetProblem(lambda, nu, bc)
  call GetWorkspace(this)

end subroutine SetTopLevelProblem_CI

!-------------------------------------------------------------------------------
!> Initialize fine-level problem with a combination of a constant isotropic
!> diffusivity and a constant isotropic spectral diffusivity

subroutine SetTopLevelProblem_CI_svv(this, lambda, nu, nu_svv, bc)
  class(PMG_Level), intent(inout) :: this
  real(RNP),        intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),        intent(in)    :: nu     !< diffusivity
  real(RNP),        intent(in)    :: nu_svv !< spectral diffusivity
  character,        intent(in)    :: bc(:)  !< boundary conditions

  call this % elliptic_op % SetProblem(lambda, nu, nu_svv, bc)
  call GetWorkspace(this)

end subroutine SetTopLevelProblem_CI_svv

!-------------------------------------------------------------------------------
!> Initialize fine-level problem with variable isotropic diffusivity

subroutine SetTopLevelProblem_VI(this, lambda, nu, bc)
  class(PMG_Level), intent(inout) :: this
  real(RNP),        intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP),        intent(in)    :: nu(0:,0:,0:,:) !< diffusivity
  character,        intent(in)    :: bc(:)          !< boundary conditions

  call this % elliptic_op % SetProblem(lambda, nu, bc)
  call GetWorkspace(this)

end subroutine SetTopLevelProblem_VI

!-------------------------------------------------------------------------------
!> Initialize coarse-level problem

subroutine SetCoarseProblem(this, fine)
  class(PMG_Level), intent(inout) :: this
  class(PMG_Level), intent(in)    :: fine !< next finer level

  real(RNP), allocatable, save :: nu_vi(:,:,:,:)

  if (allocated(fine % elliptic_op % nu_ci)) then

    call this % elliptic_op % SetProblem( fine % elliptic_op % lambda, &
                                          fine % elliptic_op % nu_ci,  &
                                          fine % elliptic_op % bc      )

  else if (allocated(fine % elliptic_op % nu_vi)) then

    !$omp single
    allocate(nu_vi(0:this%po, 0:this%po, 0:this%po, this%ne))
    !$omp end single

   !call fine % Truncate    ( fine % elliptic_op % nu_vi, nu_vi )
    call fine % Interpolate ( fine % elliptic_op % nu_vi, nu_vi )
    call this % elliptic_op % SetProblem( fine % elliptic_op % lambda, &
                                          nu_vi,                       &
                                          fine % elliptic_op % bc      )
    !$omp barrier
    !$omp master
    deallocate(nu_vi)
    !$omp end master

  end if

  call GetWorkspace(this)

end subroutine SetCoarseProblem

!-------------------------------------------------------------------------------
!> Allocate workspace

subroutine GetWorkspace(this)
  class(PMG_Level), intent(inout) :: this

  !$omp barrier
  !$omp single
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
      allocate(this % f(0:po, 0:po, 0:po, ne))
      allocate(this % v(0:po, 0:po, 0:po, ne))
    end if

  end associate
  !$omp end single

end subroutine GetWorkspace

!-------------------------------------------------------------------------------
!> Delete workspace

subroutine FreeWorkspace(this)
  class(PMG_Level), intent(inout) :: this

  !$omp barrier
  !$omp single
  if (allocated(this % u)) deallocate(this % u)
  if (allocated(this % f)) deallocate(this % f)
  if (allocated(this % v)) deallocate(this % v)
  !$omp end single

end subroutine FreeWorkspace

!===============================================================================
! Transfer operators

!-------------------------------------------------------------------------------
!> Prolongation: coarse-to-fine interpolation

subroutine Prolongate(this, uc, uf)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uc(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uf(0:,0:,0:,:) !< fine (child) mesh variable

  call TPO_AAA(this%p2f_op, uc, uf)

end subroutine Prolongate

!-------------------------------------------------------------------------------
!> Fine-to-coarse interpolation

subroutine Interpolate(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA(this%i2c_op, uf, uc)

end subroutine Interpolate

!-------------------------------------------------------------------------------
!> Fine-to-coarse restriction -- no weighting / no assembly !!

subroutine Restrict(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA(this%r2c_op, uf, uc)

end subroutine Restrict

!-------------------------------------------------------------------------------
!> Fine-to-coarse restriction -- no weighting / no assembly !!

subroutine Truncate(this, uf, uc)
  class(PMG_Level), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< coarse (parent) mesh variable

  call TPO_AAA(this%t2c_op, uf, uc)

end subroutine Truncate

!===============================================================================
! Helpers

!-------------------------------------------------------------------------------
!> Build coarse-to-fine transfer operators

subroutine Build_C2F_TransferOps(this, fine)
  class(PMG_Level), intent(inout) :: this !< current (coarse) level
  class(PMG_Level), intent(in)    :: fine !< fine level

  real(RNP) :: xc(0:this%po), xf(0:fine%po)
  integer   :: i, k, pc, pf

  !$omp single

  pc = this % po
  pf = fine % po
  xc = this % elliptic_op % eop % x
  xf = fine % elliptic_op % eop % x

  ! prolongation to fine level ...............................................

  allocate(this % p2f_op(0:pf, 0:pc))
  do k = 0, pc
  do i = 0, pf
    this % p2f_op(i,k) = LobattoPolynomial(k, xc, xf(i))
  end do
  end do

  !$omp end single

end subroutine Build_C2F_TransferOps

!-------------------------------------------------------------------------------
!> Build fine-to-coarse transfer operators

subroutine Build_F2C_TransferOps(this, pc)
  class(PMG_Level), intent(inout) :: this !< current (fine) level
  integer,          intent(in)    :: pc   !< polynomial order of coarse level

  real(RNP), allocatable :: VL(:,:), VL_inv(:,:)
  real(RNP) :: xf(0:this%po), xc(0:pc)
  integer   :: i, k, pf

  !$omp single

  pf = this % po
  xf = this % elliptic_op % eop % x
  xc = LobattoPoints(pc)

  ! interpolation to child level ...............................................

  allocate(this % i2c_op(0:pc,0:pf))
  do k = 0, pf
  do i = 0, pc
    this % i2c_op(i,k) = LobattoPolynomial(k, xf, xc(i))
  end do
  end do

  ! restriction to child level .................................................

  allocate(this % r2c_op(0:pc,0:pf))
  do k = 0, pf
  do i = 0, pc
    this % r2c_op(i,k) = LobattoPolynomial(i, xc, xf(k))
  end do
  end do

  ! truncation to child level ..................................................

  ! initialize with fine-to-coarse interpolation
  allocate(this % t2c_op(0:pc,0:pf), source = this % i2c_op)

  ! prepend order reduction, if possible
  associate(eop => this % elliptic_op % eop, t2c_op => this % t2c_op)
    if (eop % Has_Legendre_VDM()) then
      allocate(VL(0:pf,0:pf), VL_inv(0:pf,0:pf))
      call eop % Get_Legendre_VDM(VL)
      call eop % Get_Inverse_Legendre_VDM(VL_inv)
      t2c_op = matmul(t2c_op, matmul(VL(:,0:pc), VL_inv(0:pc,:)))
    end if
  end associate

  !$omp end single

end subroutine Build_F2C_TransferOps

!===============================================================================

end module CART__Elliptic_PMG_Level
