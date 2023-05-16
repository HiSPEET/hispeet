!> summary: 1D DG elliptic operator: evaluation with regular mesh and constant ν
!> author:  Joerg Stiller
!> date:    2023/03/15
!> license: Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__1D) MP_Eval_RC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Evaluation with regular (equidistant) mesh
  !>
  !> Usage
  !>   1) `f` and `bv` given:  computation of the residual, `r = f - Au`
  !>   2) `f` and `bv` absent: evaluation of the homogeneous operator, `r = Au`
  !>
  !> Boundary conditions and values
  !>   - Dirichlet:  `bc = 'D',  bv = u`
  !>   - Neumann:    `bc = 'N',  bv = q = ν∂u/∂x`
  !>
  !> Note that the flux `q` is aligned with the x-direction and not with the
  !> normal!
  !>
  !> The operator is applied only to elements for which `mask` is true.
  !> For other elements, the result is set to zero.

  module subroutine Eval_RC(this, bc, mask, dx, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_1D),   intent(in)  :: this
    character,                       intent(in)  :: bc(2)   !< BC {'D','N','P'}
    logical,                         intent(in)  :: mask(:) !< element mask
    real(RNP),                       intent(in)  :: dx      !< ∆xᵉ
    real(RNP),                       intent(in)  :: lambda  !< λ
    real(RNP),                       intent(in)  :: nu      !< ν = νᵖ+νˢ
    real(RNP), contiguous,           intent(in)  :: u(0:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r(0:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f(0:,:) !< RHS
    real(RNP),             optional, intent(in)  :: bv(2)   !< boundary values

    ! local variables ..........................................................

    real(RNP), allocatable, save :: jmp_u(:) ! [u]
    real(RNP), allocatable, save :: jmp_q(:) ! [q]
    real(RNP), allocatable, save :: avg_q(:) ! {q}

    ! 1D standard operators
    real(RNP) :: As(0:this%eop%po, 0:this%eop%po) ! diffusion Dᵀ(ν+νˢQ)D
    real(RNP) :: Bs(0:this%eop%po, 0:this%eop%po) ! "flux" (ν+νˢQ)D

    real(RNP) :: mu, nu_p, nu_s
    real(RNP) :: cl, cr, g0, g1, gh, ql, qr
    integer   :: po, ne
    integer   :: e, i
    logical   :: has_bv, has_f, hybrid

    if (nu <= 0 .or. size(r) == 0) then
      r = 0
      return
    end if

    associate(eop => this % eop, M => this % eop % w)

      ! initialization .........................................................

      has_bv   = present(bv)
      has_f    = present(f)
      hybrid   = eop % hybrid

      nu_p = this % PhysicalDiffusivity(nu)
      nu_s = this % SpectralDiffusivity(nu)

      po = eop  % po
      ne = size(u,2)
      mu = eop % PenaltyFactor(dx)
      g0 = dx / 2
      g1 = 2  / dx
      gh = g1 / (4 * mu * nu)

      ! computation of 1D standard diffusion and standard flux operator
      if (nu_s > 0) then
        call eop % Get_SVV_StandardStiffnessMatrix(As)
        call eop % Get_SVV_StandardDiffMatrix(Bs)
        As = nu_s * As
        Bs = nu_s * Bs
      else
        As = ZERO
        Bs = ZERO
      end if

      As = As + nu_p * eop%L
      Bs = Bs + nu_p * eop%D

      !$omp master
      allocate(jmp_u(0:ne), source = ZERO)
      allocate(jmp_q(0:ne), source = ZERO)
      allocate(avg_q(0:ne), source = ZERO)
      !$omp end master

      ! application of interior operator .......................................

      !$omp do
      do e = 1, ne
        if (mask(e)) then
          r(:,e) = lambda * g0 * M * u(:,e)  +  g1 * matmul(As, u(:,e))
        else
          r(:,e) = ZERO
        end if
      end do

      ! interior fluxes ........................................................

      !$omp do
      do e = 1, ne-1
        jmp_u(e) = u(po,e) - u(0,e+1)
        ql = g1 * dot_product(Bs(po,:), u(:,e  ))
        qr = g1 * dot_product(Bs( 0,:), u(:,e+1))
        jmp_q(e) = (ql - qr)
        avg_q(e) = (ql + qr) * HALF
      end do

      ! boundary fluxes ........................................................

      !$omp master

      ! left boundary
      select case(bc(1))

      case('D')
        if (has_bv) then
          jmp_u(0) = 2 * (bv(1) - u(0,1))
        else
          jmp_u(0) = 2 * (      - u(0,1))
        end if
        avg_q(0) = g1 * dot_product(Bs(0,:), u(:,1))

      case('N')
        if (has_bv) then
          avg_q(0) = bv(1)
        end if

      case('P')
        jmp_u(0) = u(po,ne) - u(0,1)
        ql = g1 * dot_product(Bs(po,:), u(:,ne))
        qr = g1 * dot_product(Bs( 0,:), u(:, 1))
        jmp_q(0) = (ql - qr)
        avg_q(0) = (ql + qr) * HALF
      end select

      ! right boundary
      select case(bc(2))

      case('D')
        if (has_bv) then
          jmp_u(ne) = 2 * (u(po,ne) - bv(2))
        else
          jmp_u(ne) = 2 * (u(po,ne)        )
        end if
        avg_q(ne) = g1 * dot_product(Bs(po,:), u(:,ne))

      case('N')
        if (has_bv) then
          avg_q(ne) = bv(2)
        end if

      case('P')
        jmp_u(ne) = jmp_u(0)
        jmp_q(ne) = jmp_q(0)
        avg_q(ne) = avg_q(0)
      end select

      ! apply fluxes and RHS ...................................................

      !$omp do
      do e = 1, ne
        if (mask(e)) then

          ! - {ν∂v/∂x}[u]
          cl = -g1/2 * jmp_u(e-1)
          cr = -g1/2 * jmp_u(e  )
          do i = 0, po
            r(i,e) = r(i,e) + cl * Bs(0,i) + cr * Bs(po,i)
          end do

          ! - [v]{ν∂u/∂x}
          r( 0,e) = r( 0,e) + avg_q(e-1)
          r(po,e) = r(po,e) - avg_q(e  )

          ! + μ⟨ν⟩[v][u]
          r( 0,e) = r( 0,e) - mu * (nu_p + nu_s) * jmp_u(e-1)
          r(po,e) = r(po,e) + mu * (nu_p + nu_s) * jmp_u(e  )

          ! - 1/4μ⟨ν⟩ [ν∂v/∂x][ν∂u/∂x]
          if (hybrid) then
            cl =  gh * jmp_q(e-1)
            cr = -gh * jmp_q(e  )
            do i = 0, po
              r(i,e) = r(i,e) + cl * Bs( 0,i) + cr * Bs(po,i)
            end do
          end if

          if (has_f) then
            r(:,e) = f(:,e) - r(:,e)
          end if

        end if
      end do

      ! clean-up ...............................................................

      !$omp master
      deallocate(jmp_u, jmp_q, avg_q)
      !$omp end master

    end associate

  end subroutine Eval_RC

  !=============================================================================

end submodule MP_Eval_RC
