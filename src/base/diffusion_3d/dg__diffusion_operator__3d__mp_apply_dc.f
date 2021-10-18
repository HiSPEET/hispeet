!> summary:  3D DG diffusion operator: application with deformed mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/10/6
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D:MP_Apply) MP_Apply_DC
  use TPO__Diffusion__3D
  use Mesh_Metrics__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

! NEXT LINE IS REQUIRED TO INCORPORATE THE EXPERIMENTAL TPO KERNEL
integer, parameter :: RWP = RNP

contains

  !-----------------------------------------------------------------------------
  !> Application with irregular (deformed) mesh and constant diffusivity

  module subroutine Apply_DC(this, u, v, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result
    real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: tr_buf

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of u and q_n

    integer :: po, ne, ng, np

    associate( mesh   => this % sem % mesh           &
             , Jd     => this % sem % metrics % Jd   &
             , G      => this % sem % metrics % G    &
             , a      => this % sem % metrics % a    &
             , Ji_n   => this % sem % metrics % Ji_n &
             , lambda => this % lambda               &
             , nu     => this % nu_pc                &
             , eop    => this % eop                  )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,2))
      tr_buf = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Diffusion_DLCI_Gen_RWP( &
                          eop%w, eop%D, Jd, G, lambda, nu, u, v &
                        , Ji_n                                  &
                        , ub = tr(:,:,:,:,1)                    &
                        , qb = tr(:,:,:,:,2)                    )

      ! transfer traces ........................................................

      call tr_buf % Transfer(mesh, tr, tag=1000)
      call tr_buf % Merge(tr)

      call ApplyBoundaryConditions( mesh, this%bc          &
                                  , tr_u  = tr(:,:,:,:,1)  &
                                  , tr_qn = tr(:,:,:,:,2)  )

      ! add fluxes .............................................................

      call AddFluxes( mesh, eop, a, Ji_n, nu &
                    , tr_u  = tr(:,:,:,:,1)  &
                    , tr_qn = tr(:,:,:,:,2)  &
                    , f     = f              &
                    , v     = v              )

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_DC

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(mesh, eop, a, Ji_n, nu, tr_u, tr_qn, f, v)

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    class(DG_ElementOperators_1D), intent(in) :: eop !< ID-DG element operators

    real(RNP), intent(in)    :: a(0:,0:,:,:)      !< area coeff @ element faces
    real(RNP), intent(in)    :: Ji_n(0:,0:,:,:,:) !< J⁻¹⋅n      @ element faces
    real(RNP), intent(in)    :: nu                !< diffusivity
    real(RNP), intent(in)    :: tr_u (0:,0:,:,:)  !< u nᵢ       @ element faces
    real(RNP), intent(in)    :: tr_qn(0:,0:,:,:)  !< q_n        @ element faces
    real(RNP), intent(inout) :: v(0:,0:,0:,:)     !< result

    real(RNP), optional, intent(in) :: f(0:,0:,0:,:) !< RHS

    ! local variables ..........................................................

    real(RNP), dimension(0:eop%po, 0:eop%po, 3) :: B
    real(RNP), dimension(0:eop%po, 0:eop%po)    :: Mf, Aq, Ju
    real(RNP) :: mu_nu, tmp
    integer   :: e, i, j, k, l, m
    logical   :: present_f

    associate( P  => eop % po, &
               Ms => eop % w,  &
               Ds => eop % D   )

      ! auxiliaries ............................................................

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      present_f = present(f)

      ! add fluxes .............................................................

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east)  . . . . . . . . . . . . . . . . . . . . . .

          do l = 1, 2

            ! Ju = n·[u], Aq = n·{ν∇u}
            call GetElementBoundaryFluxes(element, e, l, tr_u, tr_qn, Ju, Aq)

            mu_nu = eop % PenaltyFactor(mesh % dx_mean(l,e)) * nu

            i = (l-1) * P

            do k = 0, P
            do j = 0, P
              tmp = HALF * nu * a(j,k,l,e) * Ju(j,k)
              B(j,k,1) = tmp * Ji_n(j,k,l,e,1)
              B(j,k,2) = tmp * Ji_n(j,k,l,e,2)
              B(j,k,3) = tmp * Ji_n(j,k,l,e,3)
            end do
            end do

            do k = 0, P
            do j = 0, P

              do m = 0, P
                v(m,j,k,e) = v(m,j,k,e) - Mf(j,k) * Ds(i,m) * B(j,k,1)
              end do

              tmp = a(j,k,l,e) * (Aq(j,k) - mu_nu * Ju(j,k))
              do m = 0, P
                tmp = tmp + Ds(m,j) * B(m,k,2) + Ds(m,k) * B(j,m,3)
              end do
              v(i,j,k,e) = v(i,j,k,e) - Mf(j,k) * tmp

            end do
            end do

          end do

          ! faces 3+4 (south + north)  . . . . . . . . . . . . . . . . . . . . .

          do l = 3, 4

            call GetElementBoundaryFluxes(element, e, l, tr_u, tr_qn, Ju, Aq)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(l,e)) * nu
            j = (l-3) * P

            do k = 0, P
            do i = 0, P
              tmp = HALF * nu * a(i,k,l,e) * Ju(i,k)
              B(i,k,1) = tmp * Ji_n(i,k,l,e,1)
              B(i,k,2) = tmp * Ji_n(i,k,l,e,2)
              B(i,k,3) = tmp * Ji_n(i,k,l,e,3)
            end do
            end do

            do k = 0, P
            do i = 0, P

              do m = 0, P
                v(i,m,k,e) = v(i,m,k,e) - Mf(i,k) * Ds(j,m) * B(i,k,2)
              end do

              tmp = a(i,k,l,e) * (Aq(i,k) - mu_nu * Ju(i,k))
              do m = 0, P
                tmp = tmp + Ds(m,i) * B(m,k,1) + Ds(m,k) * B(i,m,3)
              end do
              v(i,j,k,e) = v(i,j,k,e) - Mf(i,k) * tmp

            end do
            end do

          end do

          ! faces 5+6 (bottom + top) . . . . . . . . . . . . . . . . . . . . . .

          do l = 5, 6

            call GetElementBoundaryFluxes(element, e, l, tr_u, tr_qn, Ju, Aq)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(l,e)) * nu
            k = (l-5) * P

            do j = 0, P
            do i = 0, P
              tmp = HALF * nu * a(i,j,l,e) * Ju(i,j)
              B(i,j,1) = tmp * Ji_n(i,j,l,e,1)
              B(i,j,2) = tmp * Ji_n(i,j,l,e,2)
              B(i,j,3) = tmp * Ji_n(i,j,l,e,3)
            end do
            end do

            do j = 0, P
            do i = 0, P

              do m = 0, P
                v(i,j,m,e) = v(i,j,m,e) - Mf(i,j) * Ds(k,m) * B(i,j,3)
              end do

              tmp = a(i,j,l,e) * (Aq(i,j) - mu_nu * Ju(i,j))
              do m = 0, P
                tmp = tmp + Ds(m,i) * B(m,j,1) + Ds(m,j) * B(i,m,2)
              end do
              v(i,j,k,e) = v(i,j,k,e) - Mf(i,j) * tmp

            end do
            end do

          end do

          ! RHS  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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

! THIS IS AN EXPERIMENTAL VERSION OF THE GENERIC TPO KERNEL WITH FLUXES

!> summary:  3d generic curvilinear element diffusion operator (DLCI)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/07/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Diffusion_DLCI_Gen_RWP(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                      Ji_n, ub, qb)

  !---------------------------------------------------------------------------
  ! Arguments

  real(RWP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
  real(RWP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
  real(RWP), intent(in)  :: Jd(:,:,:,:) !< Jacobian determinant  (np,np,np,ne)
  real(RWP), intent(in)  :: G(:,:,:,:,:)!< Laplacian metrics     (np,np,np,ne,6)
  real(RWP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RWP), intent(in)  :: nu          !< diffusivity
  real(RWP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

  real(RWP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
  real(RWP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
  real(RWP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

  !---------------------------------------------------------------------------
  ! Local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, r, s, t, z

  real(RWP) :: tmp
  integer   :: np, ne
  integer   :: e, f, i, j, k, p
  logical   :: get_traces

  !---------------------------------------------------------------------------
  ! Initialization

  np = size(Ms)
  ne = size(u,4)

  get_traces = present(Ji_n) .and. present(ub) .and. present(qb)

  ! element mass matrix

  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  !-----------------------------------------------------------------------------
  ! Evaluation

  !$omp do private(e)
  do e = 1, ne

    ! lambda M u ...............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
    end do
    end do
    end do

    ! standard derivatives of uᵉ ...............................................

    ! direction 1: r = [I x I x Dˢ] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(i,p) * u(p,j,k,e)
      end do
      r(i,j,k) = tmp
    end do
    end do
    end do

    ! direction 2: s = [I x Dˢ x I] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(j,p) * u(i,p,k,e)
      end do
      s(i,j,k) = tmp
    end do
    end do
    end do

    ! direction 3: t = [Dˢ x I x I] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(k,p) * u(i,j,p,e)
      end do
      t(i,j,k) = tmp
    end do
    end do
    end do

    ! traces on element boundary ...............................................

    if (get_traces) then

      ! ub = u,  qb = n⋅ν∇u  @ Γ₁ ∪ Γ₂

      i = 1
      do f = 1, 2
        do k = 1, np
        do j = 1, np
          ub(j,k,f,e) = u(i,j,k,e)
          qb(j,k,f,e) = nu * ( Ji_n(j,k,f,e,1) * r(i,j,k) &
                             + Ji_n(j,k,f,e,2) * s(i,j,k) &
                             + Ji_n(j,k,f,e,3) * t(i,j,k) )
        end do
        end do
        i = np
      end do

      ! ub = u,  qb = n⋅ν∇u  @ Γ₃ ∪ Γ₄

      j = 1
      do f = 3, 4
        do k = 1, np
        do i = 1, np
          ub(i,k,f,e) = u(i,j,k,e)
          qb(i,k,f,e) = nu * ( Ji_n(i,k,f,e,1) * r(i,j,k) &
                             + Ji_n(i,k,f,e,2) * s(i,j,k) &
                             + Ji_n(i,k,f,e,3) * t(i,j,k) )
        end do
        end do
        j = np
      end do

      ! ub = u,  qb = n⋅ν∇u  @ Γ₅ ∪ Γ₆

      k = 1
      do f = 5, 6
        do j = 1, np
        do i = 1, np
          ub(i,j,f,e) = u(i,j,k,e)
          qb(i,j,f,e) = nu * ( Ji_n(i,j,f,e,1) * r(i,j,k) &
                             + Ji_n(i,j,f,e,2) * s(i,j,k) &
                             + Ji_n(i,j,f,e,3) * t(i,j,k) )
        end do
        end do
        k = np
      end do

    end if

    ! completion: v += ∇ˢ⋅(ν M G ⋅ ∇ˢuᵉ) .......................................

    ! z = ν M G(1,:) ⋅ ∇ˢuᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,1) * r(i,j,k) &
                                 + G(i,j,k,e,2) * s(i,j,k) &
                                 + G(i,j,k,e,3) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [I x I x (Dˢ)ᵀ] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(p,i) * z(p,j,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + tmp
    end do
    end do
    end do

    ! z = ν M G(2,:) ⋅ ∇ˢuᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,2) * r(i,j,k) &
                                 + G(i,j,k,e,4) * s(i,j,k) &
                                 + G(i,j,k,e,5) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [I x (Dˢ)ᵀ x I] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(p,j) * z(i,p,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + tmp
    end do
    end do
    end do

    ! z = ν M G(3,:) ⋅ ∇ˢuᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,3) * r(i,j,k) &
                                 + G(i,j,k,e,5) * s(i,j,k) &
                                 + G(i,j,k,e,6) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [(Dˢ)ᵀ x I x I] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds(p,k) * z(i,j,p)
      end do
      v(i,j,k,e) = v(i,j,k,e) + tmp
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_Diffusion_DLCI_Gen_RWP

end submodule MP_Apply_DC
