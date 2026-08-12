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

!> summary: 1D DG elliptic operator: evaluation with regular mesh and variable ν
!> author:  Joerg Stiller
!> date:    2023/03/15
!===============================================================================

submodule(DG__Elliptic_Operator__1D) MP_Eval_RV
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Evaluation with regular (equidistant) mesh
  !>
  !> Usage
  !>   1) `f` and `bv` given:  computation of the residual, `r = f - Au`
  !>   2) `f` and `bv` absent: evaluation of the homogeneous operator, `r = Au`
  !>
  !> Note that `bv` needs to given for the residual even in periodic case, since
  !> it triggers also the correct treatment of interior boundaries, i.e. between
  !> active and frozen elements
  !>
  !> Boundary conditions and values
  !>   - Dirichlet:  `bc = 'D',  bv = u`
  !>   - Neumann:    `bc = 'N',  bv = ∂u/∂x`
  !>
  !> The operator is applied only to elements for which `mask` is true.
  !> For other elements, the result is set to zero.

  module subroutine Eval_RV(this, bc, mask, dx, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_1D),   intent(in)  :: this
    character,                       intent(in)  :: bc(2)    !< BC {'D','N','P'}
    logical,                         intent(in)  :: mask(:)  !< element mask
    real(RNP),                       intent(in)  :: dx       !< ∆xᵉ
    real(RNP),                       intent(in)  :: lambda   !< λ
    real(RNP), contiguous,           intent(in)  :: nu(0:,:) !< ν = νᵖ
    real(RNP), contiguous,           intent(in)  :: u (0:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r (0:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f (0:,:) !< RHS
    real(RNP),             optional, intent(in)  :: bv(2)    !< boundary values

    ! local variables ..........................................................

    real(RNP), allocatable, save :: dx_u(:,:) ! ∂uᵉ/∂x
    real(RNP), allocatable, save :: nu_max(:) ! ⟨ν⟩
    real(RNP), allocatable, save :: jmp_u(:)  ! [u]
    real(RNP), allocatable, save :: jmp_q(:)  ! [q]
    real(RNP), allocatable, save :: avg_q(:)  ! {q}

    real(RNP) :: g0, g1, gh, mu
    real(RNP) :: cl, cr, ql, qr
    integer   :: po, ne
    integer   :: e, el, er, i
    logical   :: has_bv, has_f, hybrid

    if (size(r) == 0) return

    associate( eop => this % eop     &
             , M   => this % eop % w &
             , D   => this % eop % D )

      ! initialization .........................................................

      has_bv   = present(bv)
      has_f    = present(f)
      hybrid   = eop % hybrid

      po = eop % po
      ne = size(u,2)
      mu = eop % PenaltyFactor(dx)
      g0 = dx / 2
      g1 = 2  / dx
      gh = g1 / (4 * mu)

      !$omp master
      allocate(dx_u(0:po,1:ne))
      allocate(nu_max(0:ne), source = ZERO)
      allocate(jmp_u (0:ne), source = ZERO)
      allocate(jmp_q (0:ne), source = ZERO)
      allocate(avg_q (0:ne), source = ZERO)
      !$omp end master

      ! application of interior operator .......................................

      !$omp do
      do e = 1, ne
        if (mask(e)) then
          dx_u(:,e) = g1 * matmul(D, u(:,e))
          r(:,e) = lambda * g0 * M * u(:,e) + matmul(M * nu(:,e) * dx_u(:,e), D)
        else
          dx_u( 0,e) = g1 * dot_product(D( 0,:), u(:,e))
          dx_u(po,e) = g1 * dot_product(D(po,:), u(:,e))
          r(:,e) = ZERO
        end if
      end do

      ! interior fluxes ........................................................

      !$omp do
      do e = 1, ne-1
        nu_max(e) = max(nu(po,e), nu(0,e+1))
        jmp_u (e) = u(po,e) - u(0,e+1)
        ql = nu(po,e  ) * dx_u(po,e  )
        qr = nu( 0,e+1) * dx_u( 0,e+1)
        jmp_q(e) = (ql - qr)
        avg_q(e) = (ql + qr) * HALF
      end do

      ! boundary fluxes ........................................................

      !$omp master

      ! left boundary
      select case(bc(1))

      case('D')
        nu_max(0) = nu(0,1)
        if (has_bv) then
          jmp_u(0) = 2 * (bv(1) - u(0,1))
        else
          jmp_u(0) = 2 * (      - u(0,1))
        end if
        avg_q(0) = nu(0,1) * dx_u(0,1)

      case('N')
        nu_max(0) = nu(0,1)
        if (has_bv) then
          avg_q(0) = nu(0,1) * bv(1)
        end if

      case('P')
        nu_max(0) = max(nu(po,ne), nu(0,1))
        jmp_u (0) = u(po,ne) - u(0,1)
        ql = nu(po,ne) * dx_u(po,ne)
        qr = nu( 0, 1) * dx_u( 0, 1)
        jmp_q(0) = (ql - qr)
        avg_q(0) = (ql + qr) * HALF
      end select

      ! right boundary
      select case(bc(2))

      case('D')
        nu_max(ne) = nu(po,ne)
        if (has_bv) then
          jmp_u(ne) = 2 * (u(po,ne) - bv(2))
        else
          jmp_u(ne) = 2 * (u(po,ne)        )
        end if
        avg_q(ne) = nu(po,ne) * dx_u(po,ne)

      case('N')
        nu_max(ne) = nu(po,ne)
        if (has_bv) then
          avg_q(ne) = nu(po,ne) * bv(2)
        end if

      case('P')
        nu_max(ne) = nu_max(0)
        jmp_u (ne) = jmp_u (0)
        jmp_q (ne) = jmp_q (0)
        avg_q (ne) = avg_q (0)
      end select

      !$omp end master
      !$omp barrier

      ! correct fluxes at interior boundaries ..................................

      !$omp do
      do e = 0, ne-1

        if (e == 0) then
          if (bc(1) /= 'P') cycle ! skip boundary
          el = ne
          er = 1
        else
          el = e
          er = e + 1
        end if

        if (mask(el) .and. .not. mask(er)) then
          ! left element active, right frozen
          nu_max(e) = nu(po,el)
          if (has_bv) then
            ! use frozen element data for bv
            jmp_u(e) = 2 * (u(po,el) - u(0,er))
          else
            jmp_u(e) = 2 * (u(po,el)          )
          end if
          avg_q(e) = nu(po,el) * dx_u(po,el)
          jmp_q(e) = 0

        else if (.not. mask(el) .and. mask(er)) then
          ! left element frozen, right active
          nu_max(e) = nu(0,er)
          if (has_bv) then
            ! use frozen element data for bv
            jmp_u(e) = 2 * (u(po,el) - u(0,er))
          else
            jmp_u(e) = 2 * (         - u(0,er))
          end if
          avg_q(e) = nu(0,er) * dx_u(0,er)
          jmp_q(e) = 0
        end if

      end do

      if (bc(2) == 'P') then
        nu_max(ne) = nu_max(0)
        jmp_u(ne)  = jmp_u(0)
        jmp_q(ne)  = jmp_q(0)
        avg_q(ne)  = avg_q(0)
      end if

      ! apply fluxes and RHS ...................................................

      !$omp do
      do e = 1, ne
        if (mask(e)) then

          ! - {ν∂v/∂x}[u]
          cl = -g1/2 * nu( 0,e) * jmp_u(e-1)
          cr = -g1/2 * nu(po,e) * jmp_u(e  )
          do i = 0, po
            r(i,e) = r(i,e) + cl * D(0,i) + cr * D(po,i)
          end do

          ! - [v]{ν∂u/∂x}
          r( 0,e) = r( 0,e) + avg_q(e-1)
          r(po,e) = r(po,e) - avg_q(e  )

          ! + μ⟨ν⟩[v][u]
          r( 0,e) = r( 0,e) - mu * nu_max(e-1) * jmp_u(e-1)
          r(po,e) = r(po,e) + mu * nu_max(e  ) * jmp_u(e  )

          ! - 1/4μ⟨ν⟩ [ν∂v/∂x][ν∂u/∂x]
          if (hybrid) then
            cl =  gh / nu_max(e-1) * jmp_q(e-1)
            cr = -gh / nu_max(e  ) * jmp_q(e  )
            do i = 0, po
              r(i,e) = r(i,e) + cl * D(0,i) + cr * D(po,i)
            end do
          end if

          if (has_f) then
            r(:,e) = f(:,e) - r(:,e)
          end if

        end if
      end do

      ! clean-up ...............................................................

      !$omp master
      deallocate(dx_u, nu_max, jmp_u, jmp_q, avg_q)
      !$omp end master

    end associate

  end subroutine Eval_RV

  !=============================================================================

end submodule MP_Eval_RV
