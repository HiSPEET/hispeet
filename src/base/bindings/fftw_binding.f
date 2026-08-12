!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:    Minimal interface to FFTW3
!> author:     Joerg Stiller
!> date:       2013/06/15, revised 2020/10/
!>
!> ### Note on 1D transforms
!>
!> The module provides the "Modern Fortran" interface to FFTW and defines a more
!> comfortable interface to the one-dimensional real transforms. The latter
!> supports the following discrete Fourier transforms:
!>
!> * typ = 'real' :  halfcomplex transform (FFTW_R2HC and FFTW_HC2R)
!> * typ = 'cos'  :  cosine transform (real even or DCT-I, FFTW_REDFT00)
!> * typ = 'sin'  :  sine transform (real odd or DST-I, FFTW_RODFT00)
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
!> All transforms are initialized by calling FFT_1D_InitPlan, which creates a plan
!> for the forward and backward FFTW transforms corresponding to the selected
!> type. The two plans are stored in an instance of the derived type FFT_1D_Plan.
!> For releasing the plans call Delete_FFT_1D.
!>
!> As common with FFTW, the transforms are unnormalized. To obtain a normalized
!> transform divide by 1/n in the halfcomplex (real) case, and 1/2n in the even
!> (cos) or odd (sin) case. The normalization can be applied either to the
!> forward or backward transforms.
!===============================================================================

module FFTW_Binding
  use, intrinsic :: ISO_C_Binding
  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  implicit none

  private :: RNP
  private :: Error

  !-----------------------------------------------------------------------------
  ! FFTW3 interface

  include 'fftw3.f03'

  !-----------------------------------------------------------------------------
  !> Compound holding the bidirectional plans for an FFT

  type, public :: FFT_1D
    private
    integer     :: n = 0                 !< size
    type(C_PTR) :: forward               !< plan for the forward transform
    type(C_PTR) :: backward              !< plan for the backward transform
    real(C_DOUBLE), allocatable :: ri(:) !< workspace
    real(C_DOUBLE), allocatable :: ro(:) !< workspace
  contains
    procedure :: Init_FFT_1D
    procedure :: ForwardTransform
    procedure :: BackwardTransform
    final     :: Delete_FFT_1D
  end type FFT_1D

  ! Constructor interface
  interface FFT_1D
    module procedure New_FFT_1D
  end interface

contains

!-------------------------------------------------------------------------------
!> Constructor
!>
!> The the following FFTW types are supported:
!>   * typ = 'real'  --  halfcomplex transform
!>   * typ = 'cos'   --  cosine transform
!>   * typ = 'sin'   --  sine transform

function New_FFT_1D(typ, n) result(this)
  character(len=*), intent(in) :: typ  !< FFTW type
  integer,          intent(in) :: n    !< size of the transform
  type(FFT_1D)                 :: this !< FFT_1D object

  call Init_FFT_1D(this, typ, n)

end function New_FFT_1D

!-------------------------------------------------------------------------------
!> Initialize 1D discrete Fourier transform of size n.
!>
!> The the following FFTW types are supported:
!>   * typ = 'real'  --  halfcomplex transform
!>   * typ = 'cos'   --  cosine transform
!>   * typ = 'sin'   --  sine transform

subroutine Init_FFT_1D(this, typ, n)
  class(FFT_1D),    intent(inout) :: this !< FFT_1D object
  integer,          intent(in)    :: n    !< size of the transform
  character(len=*), intent(in)    :: typ  !< FFTW type

  call Delete_FFT_1D(this)

  allocate(this % ri(n))
  allocate(this % ro(n))

  if (typ == 'real') then
    this % forward  = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_R2HC, FFTW_ESTIMATE    )
    this % backward = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_HC2R, FFTW_ESTIMATE    )
  else if (typ == 'cos') then
    this % forward  = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_REDFT00, FFTW_ESTIMATE )
    this % backward = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_REDFT00, FFTW_ESTIMATE )
  else if (typ == 'sin') then
    this % forward  = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_RODFT00, FFTW_ESTIMATE )
    this % backward = fftw_plan_r2r_1d( n, this%ri, this%ro         &
                                      , FFTW_RODFT00, FFTW_ESTIMATE )
  end if

end subroutine Init_FFT_1D

!-------------------------------------------------------------------------------
!> Delete `FFT_1D` object.

subroutine Delete_FFT_1D(this)
  type(FFT_1D) :: this

  this % n = 0

  if (C_ASSOCIATED(this % forward)) then
    call fftw_destroy_plan(this % forward)
  end if

  if (C_ASSOCIATED(this % backward)) then
    call fftw_destroy_plan(this % backward)
  end if

  if (allocated(this % ri)) deallocate(this % ri)
  if (allocated(this % ro)) deallocate(this % ro)

end subroutine Delete_FFT_1D

!-------------------------------------------------------------------------------
!> In-place 1D forward transform.

subroutine ForwardTransform(this, u)

  class(FFT_1D), intent(inout) :: this  !< FFTW object, ri and ro are changed
  real(RNP),     intent(inout) :: u(:)  !< data being transformed

  if (size(u) == this % n) then
    this % ri = u
    call fftw_execute_r2r(this % forward, this % ri, this % ro)
    u = this % ro
  else
    call Error('ForwardTransform', 'size(u) /= n', 'FFTW_Binding')
  end if

end subroutine ForwardTransform

!-------------------------------------------------------------------------------
!> In-place 1D backward transform.

subroutine BackwardTransform(this, u)

  class(FFT_1D), intent(inout) :: this  !< FFTW object, ri and ro are changed
  real(RNP),     intent(inout) :: u(:)  !< data being transformed

  if (size(u) == this % n) then
    this % ri = u
    call fftw_execute_r2r(this % forward, this % ri, this % ro)
    u = this % ro
  else
    call Error('BackwardTransform', 'size(u) /= n', 'FFTW_Binding')
  end if

end subroutine BackwardTransform

!===============================================================================

end module FFTW_Binding
