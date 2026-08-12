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

!> summary:  Set problem parameters for elliptic operator with IP/DG-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2019/02/24
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_SetProblem
  use CART__Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> Initialize problem with constant isotropic diffusivity

module subroutine SetProblem_CI(this, lambda, nu, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda !< Helmholtz parameter
  real(RNP), intent(in) :: nu     !< diffusivity
  character, intent(in) :: bc(:)  !< boundary conditions

  !$omp barrier
  !$omp master

  if (allocated(this % nu_vi )) deallocate(this % nu_vi )
  if (allocated(this % nu_hat)) deallocate(this % nu_hat)

  this % lambda = lambda
  this % nu_ci  = nu
  this % bc     = bc

  !$omp end master

  if (allocated(this % schwarz)) then
    call this % schwarz % SetProblem(this%mesh, lambda, nu, bc)
  end if

end subroutine SetProblem_CI

!-------------------------------------------------------------------------------
!> Initialize problem with a combination of constant isotropic diffusivity and
!> constant isotropic spectral diffusivity

module subroutine SetProblem_CI_svv(this, lambda, nu, nu_svv, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda !< Helmholtz parameter
  real(RNP), intent(in) :: nu     !< diffusivity
  real(RNP), intent(in) :: nu_svv !< spectral diffusivity
  character, intent(in) :: bc(:)  !< boundary conditions

  !$omp barrier
  !$omp master

  if (allocated(this % nu_vi )) deallocate(this % nu_vi )
  if (allocated(this % nu_hat)) deallocate(this % nu_hat)

  this % lambda    = lambda
  this % nu_ci     = nu
  this % nu_ci_svv = nu_svv
  this % bc        = bc

  !$omp end master

  if (allocated(this % schwarz)) then
    call this % schwarz % SetProblem(this%eop, this%mesh, lambda, nu, nu_svv, bc)
  end if

end subroutine SetProblem_CI_svv

!-------------------------------------------------------------------------------
!> Initialize problem with variable isotropic diffusivity

module subroutine SetProblem_VI(this, lambda, nu, bc)
  class(EllipticOperator3D_IP), intent(inout) :: this
  real(RNP), intent(in) :: lambda         !< Helmholtz parameter
  real(RNP), intent(in) :: nu(0:,0:,0:,:) !< diffusivity
  character, intent(in) :: bc(:)          !< boundary conditions

  ! local variables ............................................................

  type(TraceOperator), asynchronous, allocatable, save :: trace_op
  real(RNP), allocatable, save :: tr_nu(:,:,:,:)
  integer :: i, j, k, po, np

  associate(mesh => this % mesh)

    ! preparations ............................................................

    !$omp barrier
    !$omp master

    po = this % eop % po
    np = po + 1

    if (allocated(this % nu_ci))     deallocate(this % nu_ci)
    if (allocated(this % nu_ci_svv)) deallocate(this % nu_ci_svv)

    ! (re)allocate nu_vi, if necessary
    if (allocated(this % nu_vi)) then
      if (any(shape(this % nu_vi) /= shape(nu))) deallocate(this % nu_vi)
    end if
    if (.not. allocated(this % nu_vi)) then
      allocate(this % nu_vi(0:po,0:po,0:po,mesh%ne))
    end if

    ! (re)allocate nu_hat, if necessary
    if (allocated(this % nu_hat)) then
      if (any(shape(this % nu_hat) /= [np,np,mesh%nf])) then
        deallocate(this % nu_hat)
      end if
    end if
    if (.not. allocated(this % nu_hat)) then
      allocate(this % nu_hat(0:po,0:po,mesh%nf))
    end if

    ! start generating traces of nu
    allocate(tr_nu(0:po,0:po,2,mesh%nf))
    allocate(trace_op)

    !$omp end master

    ! components ...............................................................

    call trace_op % GetTrace_Start(mesh, nu, tr_nu, tag=1000)

    !$omp master
    this % lambda = lambda
    this % bc     = bc
    !$omp end master

    call SetArray(this % nu_vi, nu)

    if (allocated(this % schwarz)) then
      call this % schwarz % SetProblem(this%eop, mesh, lambda, nu, bc)
    end if

    call trace_op % GetTrace_Finish(mesh, tr_nu)

    associate(nu_hat => this % nu_hat)
      !$omp do
      do k = 1, mesh%nf
        do j = 0, po
        do i = 0, po
          nu_hat(i,j,k) = max(tr_nu(i,j,1,k), tr_nu(i,j,2,k))
        end do
        end do
      end do
    end associate

    ! clean-up .................................................................

    !$omp master
    deallocate(trace_op)
    deallocate(tr_nu)
    !$omp end master

  end associate

end subroutine SetProblem_VI

!===============================================================================

end submodule MP_SetProblem
