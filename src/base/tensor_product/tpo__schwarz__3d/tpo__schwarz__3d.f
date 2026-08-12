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

!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!> author:   Joerg Stiller
!> date:     2020/05/29
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
!===============================================================================

module TPO__Schwarz__3D
  use TPO__Schwarz__3D_I
  use TPO__Schwarz__3D_A
end module TPO__Schwarz__3D
