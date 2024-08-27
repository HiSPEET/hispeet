!> summary:  DG polynomial multigrid transfer operator
!> author:   Joerg Stiller
!> date:     2018/01/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### DG polynomial multigrid transfer operator
!===============================================================================

module CART__DG_PMG_Transfer_Operators

  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  use Standard_Element_Operators__1D
  use TPO__AAA__3D

  implicit none
  private

  public :: PMG_TransferOperators

  !-----------------------------------------------------------------------------
  !> Polynomial interpolation and restriction with discontinuous elements

  type PMG_TransferOperators
    integer :: pc = -1 !< polynomial order of coarse level
    integer :: pf = -1 !< polynomial order of fine level
    real(RNP), allocatable :: c2f_iop(:,:)  !< 1D interpolation operator
    real(RNP), allocatable :: f2c_iop(:,:)  !< 1D interpolation operator
    real(RNP), allocatable :: f2c_rop(:,:)  !< 1D restriction operator
    real(RNP), allocatable :: f2c_top(:,:)  !< 1D truncation operator
  contains
    procedure :: New         => New_PMG_TransferOperators
    procedure :: Prolongate  => C2F_Interpolation
    procedure :: Interpolate => F2C_Interpolation
    procedure :: Restrict    => F2C_Restriction
    procedure :: Truncate    => F2C_Truncation
    final     :: Delete_PMG_TransferOperators
  end type PMG_TransferOperators

contains

!-------------------------------------------------------------------------------
!> Initialization of a PMG_TransferOperators object

subroutine New_PMG_TransferOperators(this, pc, pf)
  class(PMG_TransferOperators), intent(inout) :: this
  integer, intent(in) :: pc  !< polynomial order of coarse level
  integer, intent(in) :: pf  !< polynomial order of fine level

  type(StandardElementOperators_1D) :: sop_c, sop_f
  real(RNP), allocatable :: V_c(:,:), VI_f(:,:)
  integer :: i, k

  ! clean-up ...................................................................

  call Delete_PMG_TransferOperators(this)

  ! prerequisites ..............................................................

  sop_c = StandardElementOperators_1D(pc)
  sop_f = StandardElementOperators_1D(pf)

  allocate(V_c(0:pc,0:pc), VI_f(0:pf,0:pf))
  call sop_c % Get_Legendre_VDM(V_c)
  call sop_f % Get_Inverse_Legendre_VDM(VI_f)

  this % pc = pc
  this % pf = pf

  allocate(this % c2f_iop(0:pf,0:pc))
  allocate(this % f2c_rop(0:pc,0:pf))
  allocate(this % f2c_iop(0:pc,0:pf))
  allocate(this % f2c_top(0:pc,0:pf))

  associate(xc => sop_c % x, xf => sop_f % x)

    ! coarse-to-fine interpolation .............................................

    do k = 0, pc
    do i = 0, pf
      this % c2f_iop(i,k) = LobattoPolynomial(k, xc, xf(i))
    end do
    end do

    ! fine-to-coarse interpolation .............................................

    do k = 0, pf
    do i = 0, pc
      this % f2c_iop(i,k) = LobattoPolynomial(k, xf, xc(i))
    end do
    end do

    ! fine-to-coarse restriction ...............................................

    this % f2c_rop = transpose(this % c2f_iop)

    ! fine-to-coarse truncation ................................................

    this % f2c_top = matmul(V_c, VI_f(0:pc,:))

  end associate

end subroutine New_PMG_TransferOperators

!-------------------------------------------------------------------------------
!> Finalization of a PMG_TransferOperators object

subroutine Delete_PMG_TransferOperators(this)
  type(PMG_TransferOperators), intent(inout) :: this

  if (allocated(this % c2f_iop)) deallocate(this % c2f_iop)
  if (allocated(this % f2c_rop)) deallocate(this % f2c_rop)
  if (allocated(this % f2c_iop)) deallocate(this % f2c_iop)
  if (allocated(this % f2c_top)) deallocate(this % f2c_top)

end subroutine Delete_PMG_TransferOperators

!-------------------------------------------------------------------------------
!> Prolongation: coarse-to-fine interpolation

subroutine C2F_Interpolation(this, uc, uf)
  class(PMG_TransferOperators), intent(in) :: this
  real(RNP), intent(in)  :: uc(0:,0:,0:,:) !< mesh variable of order pc
  real(RNP), intent(out) :: uf(0:,0:,0:,:) !< mesh variable of order pf

  call TPO_AAA(this%c2f_iop, uc, uf)

end subroutine C2F_Interpolation

!-------------------------------------------------------------------------------
!> Fine-to-coarse interpolation

subroutine F2C_Interpolation(this, uf, uc)
  class(PMG_TransferOperators), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable of order pf
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< mesh variable of order pc

  call TPO_AAA(this%f2c_iop, uf, uc)

end subroutine F2C_Interpolation

!-------------------------------------------------------------------------------
!> Restriction

subroutine F2C_Restriction(this, uf, uc)
  class(PMG_TransferOperators), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable of order pf
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< mesh variable of order pc

  call TPO_AAA(this%f2c_rop, uf, uc)

end subroutine F2C_Restriction

!-------------------------------------------------------------------------------
!> Truncation: spectral cut-off filtering

subroutine F2C_Truncation(this, uf, uc)
  class(PMG_TransferOperators), intent(in) :: this
  real(RNP), intent(in)  :: uf(0:,0:,0:,:) !< mesh variable of order pf
  real(RNP), intent(out) :: uc(0:,0:,0:,:) !< mesh variable of order pc

  call TPO_AAA(this%f2c_top, uf, uc)

end subroutine F2C_Truncation

!===============================================================================

end module CART__DG_PMG_Transfer_Operators
