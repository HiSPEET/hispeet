!> summary:  3D DG diffusion operator: application with regular mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D:MP_Apply) MP_Apply_RC
  use TPO__Diffusion__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application with regular (equidistant cuboidal) mesh

  module subroutine Apply_RC(this, u, v, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result
    real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: tr_buf

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of u and q_n

    ! 1D standard operators
    real(RNP) :: As(0:this%eop%po, 0:this%eop%po) ! diffusion Dᵀ(ν+νˢQ)D
    real(RNP) :: Bs(0:this%eop%po, 0:this%eop%po) ! "flux" (ν+νˢQ)D

    integer   :: po, ne, ng, np

    associate( mesh   => this % sem % mesh &
             , lambda => this % lambda     &
             , nu_p   => this % nu_pc      &
             , nu_s   => this % nu_sc      &
             , eop    => this % eop        )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

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

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,2))
      tr_buf = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Diffusion(eop%w, As, mesh%dx, lambda, ONE, u, v, Bs &
                        , ub = tr(:,:,:,:,1)                       &
                        , qb = tr(:,:,:,:,2)                       )

      ! transfer traces ........................................................

      call tr_buf % Transfer(mesh, tr, tag=1000)
      call tr_buf % Merge(tr)

      call ApplyBoundaryConditions( mesh, this%bc          &
                                  , tr_u  = tr(:,:,:,:,1)  &
                                  , tr_qn = tr(:,:,:,:,2)  )

      ! add fluxes .............................................................

      call AddFluxes( mesh, eop, Bs, nu_p, nu_s &
                    , tr_u  = tr(:,:,:,:,1)     &
                    , tr_qn = tr(:,:,:,:,2)     &
                    , f     = f                 &
                    , v     = v                 )

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_RC

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(mesh, eop, Bs, nu_p, nu_s, tr_u, tr_qn, f, v)

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    class(DG_ElementOperators_1D), intent(in) :: eop  !< ID-DG element operators

    real(RNP), intent(in) :: Bs(0:,0:)        !< 1D standard "flux" operator
    real(RNP), intent(in) :: nu_p             !< diffusivity
    real(RNP), intent(in) :: nu_s             !< spectral diffusivity
    real(RNP), intent(in) :: tr_u (0:,0:,:,:) !< u nᵢ @ element faces
    real(RNP), intent(in) :: tr_qn(0:,0:,:,:) !< q_n  @ element faces
    real(RNP), optional, intent(in) :: f(0:,0:,0:,:) !< RHS
    real(RNP), intent(inout) :: v(0:,0:,0:,:) !< result

    ! local variables ..........................................................

    real(RNP), dimension(0:eop%po, 0:eop%po) :: Mf1, Mf2, Mf3
    real(RNP), dimension(0:eop%po, 0:eop%po) :: Aq_0, Aq_P, Ju_0, Ju_P
    real(RNP), dimension(0:eop%po)           :: delta_0, delta_P
    real(RNP) :: g(3), mu(3), cp
    integer   :: i, j, k, e
    logical   :: present_f, struct

    associate( P  => eop  % po, &
               Ms => eop  % w,  &
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

      present_f = present(f)
      struct    = mesh % structured

      ! add fluxes .............................................................

      g = ONE / dx

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! direction 1
          !
          !   v = v - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₁
          !         - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₂
          !
          ! with
          !
          !   𝜑 = \ell_i(ξ) \ell_j(η) \ell_k(ζ)
          !
          ! at face 1 (ξ = -1)
          !
          !   n⋅[𝜑]         =  δ(0,i)                         =  delta_0(i)
          !   n⋅{(ν+νˢQ)∇𝜑} = -1/∆x Bs(0,i)                   = -g(1) * Bs(0,i)
          !   n⋅[u]         =  (   u⁻(j,k) -    u⁺(j,k))₁     =  Ju_0(j,k)
          !   n⋅{(ν+νˢQ)∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₁ / 2 =  Aq_0(j,k)
          !
          ! and at face 2 (ξ = +1)
          !
          !   n⋅[𝜑]         =  δ(P,i)                         =  delta_P(i)
          !   n⋅{(ν+νˢQ)∇𝜑} =  1/∆x Bs(P,i)                   =  g(1) * Bs(P,i)
          !   n⋅[u]         =  (   u⁻(j,k) -    u⁺(j,k))₂     =  Ju_P(j,k)
          !   n⋅{(ν+νˢQ)∇u} =  (n⁻⋅q⁻(j,k) - n⁺⋅q⁺(j,k))₂ / 2 =  Aq_P(j,k)

          call GetElementBoundaryFluxes(element,struct,e,1,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,2,tr_u,tr_qn,Ju_P,Aq_P)

          cp = -(nu_p + nu_s) * mu(1)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                          &
              - Mf1(j,k) * ( delta_0(i) * Aq_0(j,k)                          &
                           + delta_P(i) * Aq_P(j,k)                          &
                           + (-g(1) * Bs(0,i) + cp * delta_0(i)) * Ju_0(j,k) &
                           + ( g(1) * Bs(P,i) + cp * delta_P(i)) * Ju_P(j,k) )
          end do
          end do
          end do

          ! v = v - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₃
          !       - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₄

          call GetElementBoundaryFluxes(element,struct,e,3,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,4,tr_u,tr_qn,Ju_P,Aq_P)

          cp = -(nu_p + nu_s) * mu(2)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                          &
              - Mf2(i,k) * ( delta_0(j) * Aq_0(i,k)                          &
                           + delta_P(j) * Aq_P(i,k)                          &
                           + (-g(2) * Bs(0,j) + cp * delta_0(j)) * Ju_0(i,k) &
                           + ( g(2) * Bs(P,j) + cp * delta_P(j)) * Ju_P(i,k) )
          end do
          end do
          end do

          ! v = v - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₅
          !       - Mf ([𝜑]⋅{(ν+νˢQ)∇u} + ({(ν+νˢQ)∇𝜑} - μ(ν+νˢ)[𝜑])⋅[u])₆

          call GetElementBoundaryFluxes(element,struct,e,5,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,6,tr_u,tr_qn,Ju_P,Aq_P)

          cp = -(nu_p + nu_s) * mu(3)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                          &
              - Mf3(i,j) * ( delta_0(k) * Aq_0(i,j)                          &
                           + delta_P(k) * Aq_P(i,j)                          &
                           + (-g(3) * Bs(0,k) + cp * delta_0(k)) * Ju_0(i,j) &
                           + ( g(3) * Bs(P,k) + cp * delta_P(k)) * Ju_P(i,j) )
          end do
          end do
          end do

          if (present_f) then
            do k = 0, P
            do j = 0, P
            do i = 0, P
              v(i,j,k,e) = v(i,j,k,e) - f(i,j,k,e)
            end do
            end do
            end do
          end if

        end associate
      end do

    end associate

  end subroutine AddFluxes

  !=============================================================================

end submodule MP_Apply_RC
