!> \file       fftw_binding.f
!> \brief      Minimal interface to the real one-dimensional transforms of FFTW.
!> \author     Joerg Stiller
!> \date       2013/06/15
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!> The module provides the "Modern Fortran" interface to FFTW and defines a more
!> comfortable interface to the one-dimensional real transforms.
!> The latter supports the following discrete Fourier transforms:
!>
!> * type = 'real' :  halfcomplex transform (FFTW_R2HC and FFTW_HC2R)
!> * type = 'cos'  :  cosine transform (real even or DCT-I, FFTW_REDFT00)
!> * type = 'sin'  :  sine transform (real odd or DST-I, FFTW_RODFT00)
!>
!> The forward halfcomplex transform of a real vector u(0:n-1) is defined as
!>
!>     c(k) = n−1 = sum_{j=0}^{n-1} u(j) exp(-I2πjk/n),  k = 0, n-1
!>
!> where I is the imaginary unit and c(k) the complex Fourier coefficients.
!> As the input data is real, the coefficients possess the Hermitian symmetry,
!> i.e., c(n-k) is the conjugate of c(k). Exploiting this property the transform
!> returns the purely real vector
!>
!>     v(0:n-1) = [r(0),r(1),r(2),...,r(n/2),i((n+1)/2−1),...,i(2),i(1)]
!>
!> where r(k) and i(k) represent the real and imaginary parts of c(k),
!> respectively. Division by 2 is rounded down.
!> The backward transform reads
!>
!>     u(j) = n−1 = sum_{k=0}^{n-1} c(k) exp(I2πjk/n),  j = 0, n-1
!>
!> The discrete cosine transform DCT-1 is defined as
!>
!>     v(k) = u(0) + (−1)^k u(n−1) + 2 sum_{j=1}^{n-2} u(j) cos[πjk/(n−1)]
!>
!> for k = 0, n-1.
!>
!> The discrete sine transform DST-1 is defined as
!>
!>     v(k) = 2 sum_{j=0}^{n-1} u(j) sin[π(j+1)(k+1)/(n+1)],  k = 0, n-1.
!>
!> Note that both the DCT-1 and the DST-1 are symmetric.
!>
!> All transforms are initialized by calling FFT_InitPlan, which creates a plan
!> for the forward and backward FFTW transforms corresponding to the selected
!> type. The two plans are stored in an instance of the derived type FFT_Plan.
!> For releasing the plans call FFT_DestroyPlan.
!>
!> As common with FFTW, the transforms are unnormalized. To obtain a normalized
!> transform divide by 1/n in the halfcomplex (real) case, and 1/2n in the even
!> (cos) or odd (sin) case. The normalization can be applied either to the
!> forward or backward transforms.
!>
!> \todo
!>   * Include examples or usage instructions into description
!>   * Testing
!===============================================================================

module FFTW_Binding
  use, intrinsic :: ISO_C_Binding
  use Kind_Parameters, only: RNP
  implicit none
  public

  !-----------------------------------------------------------------------------
  !> \brief   Compound holding the bidirectional plans for an FFT
  !> \author  Joerg Stiller

  type, public :: FFT_Plan
    type(C_PTR) :: forward   !< plan for the forward transform
    type(C_PTR) :: backward  !< plan for the backward transform
  end type FFT_Plan

  !-----------------------------------------------------------------------------
  ! public procedures

  public :: FFT_InitPlan
  public :: FFT_DestroyPlan
  public :: FFT_Execute

  !-----------------------------------------------------------------------------
  ! FFTW constants

  include 'fftw3.f03'

contains

!-------------------------------------------------------------------------------
!> \brief   Initialize discrete Fourier transform of size n.
!> \author  Joerg Stiller
!>
!> \details
!> The the following FFTW types are supported:
!>   * type = 'real'  --  halfcomplex transform
!>   * type = 'cos'   --  cosine transform
!>   * type = 'sin'   --  sine transform
!>

subroutine FFT_InitPlan(n, type, fft)
  integer,          intent(in)  :: n     !< size of the transform
  character(len=*), intent(in)  :: type  !< FFTW type
  type(FFT_Plan),   intent(out) :: fft   !< FFTW plan

  real(C_DOUBLE) :: ri(n), ro(n)

  if (type == 'real') then
    fft%forward  = fftw_plan_r2r_1d(n, ri, ro, FFTW_R2HC, FFTW_ESTIMATE)
    fft%backward = fftw_plan_r2r_1d(n, ri, ro, FFTW_HC2R, FFTW_ESTIMATE)
  else if (type == 'cos') then
    fft%forward  = fftw_plan_r2r_1d(n, ri, ro, FFTW_REDFT00, FFTW_ESTIMATE)
    fft%backward = fftw_plan_r2r_1d(n, ri, ro, FFTW_REDFT00, FFTW_ESTIMATE)
  else if (type == 'sin') then
    fft%forward  = fftw_plan_r2r_1d(n, ri, ro, FFTW_RODFT00, FFTW_ESTIMATE)
    fft%backward = fftw_plan_r2r_1d(n, ri, ro, FFTW_RODFT00, FFTW_ESTIMATE)
  end if

end subroutine FFT_InitPlan

!-------------------------------------------------------------------------------
!> \brief   Destroy FFT plan.
!> \author  Joerg Stiller

subroutine FFT_DestroyPlan(fft)
  type(FFT_Plan), intent(inout) :: fft   !< FFTW plan

  if (C_ASSOCIATED(fft%forward)) then
    call fftw_destroy_plan(fft%forward)
  end if

  if (C_ASSOCIATED(fft%backward)) then
    call fftw_destroy_plan(fft%backward)
  end if

end subroutine FFT_DestroyPlan

!-------------------------------------------------------------------------------
!> \brief   Perform discrete Fourier transform.
!> \author  Joerg Stiller

subroutine FFT_Execute(plan, u)

  type(C_PTR), intent(in)    :: plan  !< FFTW plan
  real(RNP),   intent(inout) :: u(:)  !< data being transformed

  real(C_DOUBLE), dimension(size(u)) :: ri, ro

  ri = u
  call fftw_execute_r2r(plan, ri, ro)
  u = ro

end subroutine FFT_Execute

!===============================================================================

end module FFTW_Binding