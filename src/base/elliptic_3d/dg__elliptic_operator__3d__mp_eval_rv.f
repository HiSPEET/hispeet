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

!> summary:  3D DG elliptic operator: evaluation with regular mesh and
!>           variable diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/09
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Eval_RV
  use TPO__Elliptic__3D
  use Mesh__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Evaluation with regular (equidistant cuboidal) mesh
  !>
  !> Usage
  !>   1) `f` and `bv` given:  computation of the residual, `r = f - Au`
  !>   2) `f` and `bv` absent: evaluation of the homogeneous operator, `r = Au`

  module subroutine Eval_RV(this, bc, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    character,                       intent(in)  :: bc(:)       !< {P,D,N}
    real(RNP),                       intent(in)  :: lambda      !< λ
    real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous,           intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r (:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f (:,:,:,:) !< RHS
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of ν, u and q_n

    logical :: hom_bc
    integer :: po, ne, ng, np

    if (this % sem % mesh % part < 0) return

    associate( eop  => this % eop        &
             , sem  => this % sem        &
             , mesh => this % sem % mesh )

      ! initialization .........................................................

      hom_bc = .not. present(bv)

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,3))
      buf_tr = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Elliptic( eop, sem            & ! As r may contain ghost entries
                       , lambda, nu, u       & ! it is explicitly restricted to
                       , r(:,:,:,:ne)        & ! local elements since no proper
                       , nub = tr(:,:,:,:,1) & ! bound checking is performed in
                       , ub  = tr(:,:,:,:,2) & ! TPO_Elliptic.
                       , qb  = tr(:,:,:,:,3) )

      ! transfer traces and apply boundary conditions ..........................

      call buf_tr % Transfer(mesh, tr, tag=1000)

      call EnforceBoundaryConditions( this, bc, bv          &
                                    , jmp_u = tr(:,:,:,:,2) &
                                    , avg_q = tr(:,:,:,:,3) )

      call buf_tr % Merge(tr)

      ! add fluxes .............................................................

      call AddFluxes(this, hom_bc, nu, tr, r, f)

      ! frozen and ghost elements are set to zero ..............................

      call SetArray(r(:,:,:,mesh%n_elem_active+1:), ZERO)

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, buf_tr)
      !$omp end master

    end associate

  end subroutine Eval_RV

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(this, hom_bc, nu, tr, r, f)

    ! arguments ................................................................

    class(DG_EllipticOperator_3D), intent(in) :: this

    logical, intent(in) :: hom_bc !< F/T for in/homogeneous BC
    real(RNP), contiguous, intent(in) :: nu(0:,0:,0:,:)   !< diffusivity ν
    real(RNP), contiguous, intent(in) :: tr(0:,0:,:,:,:)  !< traces of ν, u, q_n
    real(RNP), contiguous, intent(inout) :: r(0:,0:,0:,:) !< result
    real(RNP), contiguous, intent(in), optional :: f(0:,0:,0:,:) !< RHS

    ! local variables ..........................................................

    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: Mf1, Mf2, Mf3
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: jmp_u_0, jmp_u_P
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: avg_q_0, avg_q_P
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: nu_max_0, nu_max_P
    real(RNP), dimension(0:this%eop%po)                :: delta_0, delta_P
    real(RNP) :: g(3), mu(3)
    real(RNP) :: cd_0, cd_P, cp_0, cp_P
    integer   :: i, j, k, e
    logical   :: residual, struct

    associate( eop  => this % eop             &
             , P    => this % eop % po        &
             , Ms   => this % eop % w         &
             , Ds   => this % eop % D         &
             , mesh => this % sem % mesh      &
             , dx   => this % sem % mesh % dx )

      ! auxiliaries ............................................................

      ! face mass matrices

      g(1) = dx(2) * dx(3) / 4
      g(2) = dx(3) * dx(1) / 4
      g(3) = dx(1) * dx(2) / 4

      do j = 0, P
      do i = 0, P
        Mf1(i,j) = g(1) * Ms(i) * Ms(j)
        Mf2(i,j) = g(2) * Ms(i) * Ms(j)
        Mf3(i,j) = g(3) * Ms(i) * Ms(j)
      end do
      end do

      ! delta function
      delta_0 = [ ONE, (ZERO, i=1,P) ]
      delta_P = [ (ZERO, i=1,P), ONE ]

      ! penalties
      mu(1) = eop % PenaltyFactor(dx(1))
      mu(2) = eop % PenaltyFactor(dx(2))
      mu(3) = eop % PenaltyFactor(dx(3))

      residual = present(f)
      struct   = mesh % structured

      ! add fluxes .............................................................

      g = ONE / dx

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! direction 1
          !
          !   r = r - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₁
          !         - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₂
          !
          ! with
          !
          !   𝜑 = \ell_i(ξ) \ell_j(η) \ell_k(ζ)
          !
          ! at face 1 (ξ = -1)
          !
          !   n⋅[𝜑]   =  δ(0,i)       =  delta_0(i)
          !   n⋅{ν∇𝜑} = -ν/∆x D(0,i)  = -g(1) * nu(0,j,k,e) * Ds(0,i)
          !
          !   n⋅[u]   =  (   u⁻(j,k) -    u⁺(j,k))₁     =  jmp_u_0(j,k)
          !   n⋅{ν∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₁ / 2 =  avg_q_0(j,k)
          !
          ! and at face 2 (ξ = +1)
          !
          !   n⋅[𝜑]   =  δ(P,i)        =  delta_P(i)
          !   n⋅{ν∇𝜑} =  ν/∆x D(P,i)   =  g(1) * nu(P,j,k,e) * Ds(P,i)
          !
          !   n⋅[u]   =  (   u⁻(j,k) -    u⁺(j,k))₂     =  jmp_u_P(j,k)
          !   n⋅{ν∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₂ / 2 =  avg_q_P(j,k)

          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 1  &
                                              , tr, nu_max_0, jmp_u_0, avg_q_0 )
          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 2  &
                                              , tr, nu_max_P, jmp_u_P, avg_q_P )

          do k = 0, P
          do j = 0, P

            cd_0 = -nu(0,j,k,e) * g(1)
            cd_P =  nu(P,j,k,e) * g(1)
            cp_0 = -nu_max_0(j,k) * mu(1)
            cp_P = -nu_max_P(j,k) * mu(1)

            do i = 0, P

              r(i,j,k,e) = r(i,j,k,e)                                              &
                - Mf1(j,k) * ( delta_0(i) * avg_q_0(j,k)                           &
                             + delta_P(i) * avg_q_P(j,k)                           &
                             + (cd_0 * Ds(0,i) + cp_0 * delta_0(i)) * jmp_u_0(j,k) &
                             + (cd_P * Ds(P,i) + cp_P * delta_P(i)) * jmp_u_P(j,k) )
            end do
          end do
          end do

          ! r = r - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₃
          !       - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₄

          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 3  &
                                              , tr, nu_max_0, jmp_u_0, avg_q_0 )
          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 4  &
                                              , tr, nu_max_P, jmp_u_P, avg_q_P )

          do k = 0, P
          do i = 0, P

            cd_0 = -nu(i,0,k,e) * g(2)
            cd_P =  nu(i,P,k,e) * g(2)
            cp_0 = -nu_max_0(i,k) * mu(2)
            cp_P = -nu_max_P(i,k) * mu(2)

            do j = 0, P

              r(i,j,k,e) = r(i,j,k,e)                                              &
                - Mf2(i,k) * ( delta_0(j) * avg_q_0(i,k)                           &
                             + delta_P(j) * avg_q_P(i,k)                           &
                             + (cd_0 * Ds(0,j) + cp_0 * delta_0(j)) * jmp_u_0(i,k) &
                             + (cd_P * Ds(P,j) + cp_P * delta_P(j)) * jmp_u_P(i,k) )
            end do
          end do
          end do

          ! r = r - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₅
          !       - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₆

          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 5  &
                                              , tr, nu_max_0, jmp_u_0, avg_q_0 )
          call this % GetElementBoundaryFluxes( element, struct, hom_bc, e, 6  &
                                              , tr, nu_max_P, jmp_u_P, avg_q_P )

          do j = 0, P
          do i = 0, P

            cd_0 = -nu(i,j,0,e) * g(3)
            cd_P =  nu(i,j,P,e) * g(3)
            cp_0 = -nu_max_0(i,j) * mu(3)
            cp_P = -nu_max_P(i,j) * mu(3)

            do k = 0, P

              r(i,j,k,e) = r(i,j,k,e)                                              &
                - Mf3(i,j) * ( delta_0(k) * avg_q_0(i,j)                           &
                             + delta_P(k) * avg_q_P(i,j)                           &
                             + (cd_0 * Ds(0,k) + cp_0 * delta_0(k)) * jmp_u_0(i,j) &
                             + (cd_P * Ds(P,k) + cp_P * delta_P(k)) * jmp_u_P(i,j) )
            end do
          end do
          end do

          if (residual) then
            ! r = f - Au
            do k = 0, P
            do j = 0, P
            do i = 0, P
              r(i,j,k,e) = f(i,j,k,e) - r(i,j,k,e)
            end do
            end do
            end do
          end if

        end associate
      end do

    end associate

  end subroutine AddFluxes

  !=============================================================================

end submodule MP_Eval_RV
