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

      ! start generation of traces .............................................

      ! provide local traces of solution u and normal diffusive flux qn
      call GetLocalTraces( np, ne, Bs, mesh%dx, u &
                         , tr_u  = tr(:,:,:,:,1)  &
                         , tr_qn = tr(:,:,:,:,2)  )

      ! start transfer of fluxes to/from remote neighbors
      call tr_buf % Transfer(mesh, tr, tag=1000)

      ! apply element diffusion operator .......................................

      call TPO_Diffusion(eop%w, As, mesh%dx, lambda, ONE, u, v)

      ! finish generation of traces ............................................

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
  !> Elementwise computation of diffusive fluxes parallel to face normals
  !>
  !> Computes the normal components of (ν+νˢQ) ∇u for all element boundary
  !> points. Entries corresponding to interior points or tangential components
  !> are set to zero.

  subroutine GetLocalTraces(np, ne, Bs, dx, u, tr_u, tr_qn)
    integer,   intent(in)    :: np             !< num points per direction
    integer,   intent(in)    :: ne             !< num elements
    real(RNP), intent(in)    :: Bs(np,np)      !< 1D standard "flux" operator
    real(RNP), intent(in)    :: dx(3)          !< element extensions
    real(RNP), intent(in)    :: u(np,np,np,ne) !< 3D scalar field
    real(RNP), intent(inout) :: tr_u (:,:,:,:) !< solution traces
    real(RNP), intent(inout) :: tr_qn(:,:,:,:) !< normal flux traces

    real(RNP) :: Bs_f(2,np)
    real(RNP) :: g(3), tmp1, tmp2
    integer   :: e, i, j, k, m

    ! initialization ...........................................................

    ! transposed normal differentiation operators for first and last point
    Bs_f(1,:) = -Bs( 1,:)  ! first
    Bs_f(2,:) =  Bs(np,:)  ! last

    ! metric coefficients
    g = 2 / dx

    !$acc data present(u,q) copyin(Bs_f,g)
    !$acc parallel
    !$acc loop gang worker

    !$omp do private(e)
    do e = 1, ne

      ! face 1+2: qn = ∓(ν+νˢQ) ∂u/∂x1 .........................................

      !$acc loop collapse(2) vector
      do k = 1, np
      do j = 1, np

        tmp1 = 0
        tmp2 = 0
        do m = 1, np
          tmp1 = tmp1 + Bs_f(1,m) * u(m,j,k,e)
          tmp2 = tmp2 + Bs_f(2,m) * u(m,j,k,e)
        end do
        tr_u (j,k,1,e) = u( 1,j,k,e)
        tr_u (j,k,2,e) = u(np,j,k,e)
        tr_qn(j,k,1,e) = g(1) * tmp1
        tr_qn(j,k,2,e) = g(1) * tmp2

      end do
      end do

      ! face 3+4: qn = ∓(ν+νˢQ) ∂u/∂x2 .........................................

      !$acc loop collapse(2) vector
      do k = 1, np
      do i = 1, np

        tmp1 = 0
        tmp2 = 0
        do m = 1, np
          tmp1 = tmp1 + Bs_f(1,m) * u(i,m,k,e)
          tmp2 = tmp2 + Bs_f(2,m) * u(i,m,k,e)
        end do
        tr_u (i,k,3,e) = u(i, 1,k,e)
        tr_u (i,k,4,e) = u(i,np,k,e)
        tr_qn(i,k,3,e) = g(2) * tmp1
        tr_qn(i,k,4,e) = g(2) * tmp2

      end do
      end do

      ! face 5+6: qn = ∓(ν+νˢQ) ∂u/∂x3 .........................................

      !$acc loop collapse(2) vector
      do j = 1, np
      do i = 1, np

        tmp1 = 0
        tmp2 = 0
        do m = 1, np
          tmp1 = tmp1 + Bs_f(1,m) * u(i,j,m,e)
          tmp2 = tmp2 + Bs_f(2,m) * u(i,j,m,e)
        end do
        tr_u (i,j,5,e) = u(i,j, 1,e)
        tr_u (i,j,6,e) = u(i,j,np,e)
        tr_qn(i,j,5,e) = g(3) * tmp1
        tr_qn(i,j,6,e) = g(3) * tmp2

      end do
      end do

    end do

    !$acc end parallel
    !$acc end data

  end subroutine GetLocalTraces

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
    real(RNP) :: g(3), mu(3)
    real(RNP) :: cd_0, cd_P, cp_0, cp_P
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

          ! direction 1: v = v + Mf₁ (-[ϕ]₁{(ν+νˢQ)∇u}₁ - {(ν+νˢQ)∇ϕ}₁[u]₁
          !                           + μ⟨ν+νˢ⟩[ϕ]₁[u]₁)
          !
          ! where, with ϕ = ℓ_ijk, at  ξ = -1
          !
          !   [ϕ]₁         = -delta_0(i)
          !   {(ν+νˢQ)∇ϕ}₁ = g(1) * Bs(0,i)
          !   [u]₁         = jmp( u (j,k,f₁) )
          !   {(ν+νˢQ)∇u}₁ = avg( q₁(j,k,f₁) )
          !
          ! and, at  ξ = +1
          !
          !   [ϕ]₁         = delta_P(i)
          !   {(ν+νˢQ)∇ϕ}₁ = g(1) * Bs(P,i)
          !   [u]₁         = jmp( u (j,k,f₂) )
          !   {(ν+νˢQ)∇u}₁ = avg( q₁(j,k,f₂) )

          call GetElementBoundaryFluxes(element,struct,e,1,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,2,tr_u,tr_qn,Ju_P,Aq_P)

          cd_0 = -g(1)
          cd_P = -g(1)
          cp_0 = -(nu_p + nu_s) * mu(1)
          cp_P = -(nu_p + nu_s) * mu(1)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                           &
              - Mf1(j,k) * ( delta_0(i) * Aq_0(j,k)                           &
                           + delta_P(i) * Aq_P(j,k)                           &
                           + (cd_0 * Bs(0,i) + cp_0 * delta_0(i)) * Ju_0(j,k) &
                           + (cd_P * Bs(P,i) + cp_P * delta_P(i)) * Ju_P(j,k) )
          end do
          end do
          end do

          ! direction 2: v = v + Mf₂ (-[ϕ]₂{(ν+νˢQ)∇u}₂ - {(ν+νˢQ)∇ϕ}₂[u]₂
          !                           + μ⟨ν+νˢ⟩[ϕ]₂[u]₂)

          call GetElementBoundaryFluxes(element,struct,e,3,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,4,tr_u,tr_qn,Ju_P,Aq_P)

          cd_0 = -g(2)
          cd_P = -g(2)
          cp_0 = -(nu_p + nu_s) * mu(2)
          cp_P = -(nu_p + nu_s) * mu(2)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                           &
              - Mf2(i,k) * ( delta_0(j) * Aq_0(i,k)                           &
                           + delta_P(j) * Aq_P(i,k)                           &
                           + (cd_0 * Bs(0,j) + cp_0 * delta_0(j)) * Ju_0(i,k) &
                           + (cd_P * Bs(P,j) + cp_P * delta_P(j)) * Ju_P(i,k) )
          end do
          end do
          end do

          ! direction 3: v = v + Mf₃ (-[ϕ]₃{(ν+νˢQ)∇u}₃ - {(ν+νˢQ)∇ϕ}₃[u]₃
          !                           + μ⟨ν+νˢ⟩[ϕ]₃[u]₃)

          call GetElementBoundaryFluxes(element,struct,e,5,tr_u,tr_qn,Ju_0,Aq_0)
          call GetElementBoundaryFluxes(element,struct,e,6,tr_u,tr_qn,Ju_P,Aq_P)

          cd_0 = -g(3)
          cd_P = -g(3)
          cp_0 = -(nu_p + nu_s) * mu(3)
          cp_P = -(nu_p + nu_s) * mu(3)

          do k = 0, P
          do j = 0, P
          do i = 0, P

            v(i,j,k,e) = v(i,j,k,e)                                           &
              - Mf3(i,j) * ( delta_0(k) * Aq_0(i,j)                           &
                           + delta_P(k) * Aq_P(i,j)                           &
                           + (cd_0 * Bs(0,k) + cp_0 * delta_0(k)) * Ju_0(i,j) &
                           + (cd_P * Bs(P,k) + cp_P * delta_P(k)) * Ju_P(i,j) )
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
