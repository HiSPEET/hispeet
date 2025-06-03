!> summary:  3D DG elliptic operator: evaluation with regular mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Eval_RC
  use TPO__Elliptic__3D_RLCI
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

  module subroutine Eval_RC(this, bc, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    character,                       intent(in)  :: bc(:)      !< {P,D,N}
    real(RNP),                       intent(in)  :: lambda     !< λ
    real(RNP),                       intent(in)  :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of u and q_n

    logical :: hom_bc
    integer :: po, ne, ng, np

    if (this % sem % mesh % part < 0) return

    associate( eop  => this % eop        &
             , mesh => this % sem % mesh )

      ! initialization .........................................................

      hom_bc = .not. present(bv)

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,2))
      buf_tr = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      ! r is restricted to local elements, as it may contain ghost
      ! entries and no bound checking is performed in TPO_Elliptic

      call TPO_Elliptic_RLCI( eop%w, eop%L, mesh%dx  &
                            , lambda, nu, u          &
                            , r(:,:,:,:ne), eop%D    &
                            , ub = tr(:,:,:,:,1)     &
                            , qb = tr(:,:,:,:,2)     )

      ! transfer traces and apply boundary conditions ..........................

      call buf_tr % Transfer(mesh, tr, tag=1000)

      call EnforceBoundaryConditions( this, bc, bv          &
                                    , jmp_u = tr(:,:,:,:,1) &
                                    , avg_q = tr(:,:,:,:,2) )

      call buf_tr % Merge(tr)

      ! add fluxes .............................................................

      call AddFluxes(mesh, eop, hom_bc, nu, tr, r, f)

      ! frozen and ghost elements are set to zero ..............................

      call SetArray(r(:,:,:,mesh%n_elem_active+1:), ZERO)

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, buf_tr)
      !$omp end master

    end associate

  end subroutine Eval_RC

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(mesh, eop, hom_bc, nu, tr, r, f)

    ! arguments ................................................................

    class(Mesh_3D),                intent(in) :: mesh !< mesh partition
    class(DG_ElementOperators_1D), intent(in) :: eop  !< ID-DG element operators

    logical,   intent(in) :: hom_bc  !< F/T for in/homogeneous BC
    real(RNP), intent(in) :: nu      !< diffusivity
    real(RNP), contiguous, intent(in)    :: tr(0:,0:,:,:,:)  !< traces of u, q_n
    real(RNP), contiguous, intent(inout) :: r(0:,0:,0:,:)    !< result
    real(RNP), contiguous, intent(in), optional :: f(0:,0:,0:,:) !< RHS

    ! local variables ..........................................................

    real(RNP), dimension(0:eop%po, 0:eop%po) :: Mf1, Mf2, Mf3
    real(RNP), dimension(0:eop%po, 0:eop%po) :: jmp_u_0, avg_q_0
    real(RNP), dimension(0:eop%po, 0:eop%po) :: jmp_u_P, avg_q_P
    real(RNP), dimension(0:eop%po)           :: delta_0, delta_P
    real(RNP) :: g(3), mu(3), cp
    integer   :: i, j, k, e
    logical   :: residual, struct

    associate( P  => eop  % po, &
               Ms => eop  % w,  &
               Ds => eop  % D,  &
               dx => mesh % dx  )

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
      do e = 1, mesh % n_elem_active
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
          !   n⋅[𝜑]   =  δ(0,i)                         =  delta_0(i)
          !   n⋅{ν∇𝜑} = -ν/∆x D(0,i)                    = -g(1) * ν Ds(0,i)
          !   n⋅[u]   =  (   u⁻(j,k) -    u⁺(j,k))₁     =  jmp_u_0(j,k)
          !   n⋅{ν∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₁ / 2 =  avg_q_0(j,k)
          !
          ! and at face 2 (ξ = +1)
          !
          !   n⋅[𝜑]   =  δ(P,i)                         =  delta_P(i)
          !   n⋅{ν∇𝜑} =  ν/∆x D(P,i)                    =  g(1) * ν Ds(P,i)
          !   n⋅[u]   =  (   u⁻(j,k) -    u⁺(j,k))₂     =  jmp_u_P(j,k)
          !   n⋅{ν∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₂ / 2 =  avg_q_P(j,k)

          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 1 &
                                       , tr, jmp_u_0, avg_q_0          )
          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 2 &
                                       , tr, jmp_u_P, avg_q_P          )

          cp = -nu * mu(1)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            r(i,j,k,e) = r(i,j,k,e)                                             &
                       - Mf1(j,k)                                               &
                          * ( delta_0(i) * avg_q_0(j,k)                         &
                            + delta_P(i) * avg_q_P(j,k)                         &
                            + (-g(1)*nu*Ds(0,i) + cp*delta_0(i)) * jmp_u_0(j,k) &
                            + ( g(1)*nu*Ds(P,i) + cp*delta_P(i)) * jmp_u_P(j,k) )
          end do
          end do
          end do

          ! r = r - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₃
          !       - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₄

          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 3 &
                                       , tr, jmp_u_0, avg_q_0          )
          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 4 &
                                       , tr, jmp_u_P, avg_q_P          )

          cp = -nu * mu(2)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            r(i,j,k,e) = r(i,j,k,e)                                             &
                       - Mf2(i,k)                                               &
                          * ( delta_0(j) * avg_q_0(i,k)                         &
                            + delta_P(j) * avg_q_P(i,k)                         &
                            + (-g(2)*nu*Ds(0,j) + cp*delta_0(j)) * jmp_u_0(i,k) &
                            + ( g(2)*nu*Ds(P,j) + cp*delta_P(j)) * jmp_u_P(i,k) )
          end do
          end do
          end do

          ! r = r - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₅
          !       - Mf ([𝜑]⋅{ν∇u} + ({ν∇𝜑} - μν[𝜑])⋅[u])₆

          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 5 &
                                       , tr, jmp_u_0, avg_q_0          )
          call GetElementBoundaryFluxes( element, struct, hom_bc, e, 6 &
                                       , tr, jmp_u_P, avg_q_P          )

          cp = -nu * mu(3)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            r(i,j,k,e) = r(i,j,k,e)                                             &
                       - Mf3(i,j)                                               &
                          * ( delta_0(k) * avg_q_0(i,j)                         &
                            + delta_P(k) * avg_q_P(i,j)                         &
                            + (-g(3)*nu*Ds(0,k) + cp*delta_0(k)) * jmp_u_0(i,j) &
                            + ( g(3)*nu*Ds(P,k) + cp*delta_P(k)) * jmp_u_P(i,j) )
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

end submodule MP_Eval_RC
