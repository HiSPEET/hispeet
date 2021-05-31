!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> The procedures evaluate a spectral tensor-product operator of the form
!>
!>     (W3xW2xW1) (S3xS2xS1) D⁻¹ (S3xS2xS1)^T
!>
!> where
!>
!>   * `S1,S2,S3` are the column matrices of the eigenvectors to the
!>      1D suboperators
!>
!>   * `D` is a diagonal matrix composed of the 1D eigenvalues
!>
!>         D = lambda  g0 I3 x I2 x I1
!>           + nu(e) ( g1 I3 x I2 x V1
!>                   + g2 I3 x V2 x I1
!>                   + g3 V3 x I2 x I1 )
!>
!>   * `W1,W2,W3` are the diagonal 1D weighting matrices
!>
!>   * metric coefficients
!>
!>         g0 = dx(1) * dx(2) * dx(3)
!>         g1 = dx(2) * dx(3) / dx(1)
!>         g2 = dx(3) * dx(1) / dx(2)
!>         g3 = dx(1) * dx(2) / dx(3)
!>
!>
!> To accomodate various types of boundary conditions, every operator
!> is given for each possible 1D boundary configuration. The boundary
!> configurations can differ from element, as specified in the parameter
!> `cfg`.
!>
!> When compiled with OpenACC and run on a GPU,
!>
!>   * cfg, D_inv, f and u must be present on the device, and
!>   * any call to this procedure must be followed by an "acc wait"
!>     before accessing f or u from CPU or another than the default
!>     accelerator queue.
!===============================================================================

module TPO__Schwarz__3D
  use TPO__Schwarz__3D_CI
  use TPO__Schwarz__3D_CA
  implicit none
  private

  public :: TPO_Schwarz

  interface TPO_Schwarz
    module procedure TPO_Schwarz_CI ! (W S  x W S  x W S ) D⁻¹ (S  x S  x S )ᵀ
    module procedure TPO_Schwarz_CA ! (W₃S₃ x W₂S₂ x W₁S₁) D⁻¹ (S₃ x S₂ x S₁)ᵀ
  end interface

end module TPO__Schwarz__3D
