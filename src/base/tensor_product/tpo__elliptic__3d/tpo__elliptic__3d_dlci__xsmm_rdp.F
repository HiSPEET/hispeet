!> summary:  3D element elliptic operator using LIBXSMM (DLCI) for real(RDP)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic__3D_DLCI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Elliptic_DLCI_XSMM

  interface TPO_Elliptic_DLCI_XSMM
    module procedure TPO_Elliptic_DLCI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element elliptic operator based on LIBXSMM (RLCI)

  subroutine TPO_Elliptic_DLCI_XSMM_RDP(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                        Ji_n, ub, qb)

    real(RDP), intent(in)  :: Ms(:)        !< 1D standard mass matrix       (np)
    real(RDP), intent(in)  :: Ds(:,:)      !< 1D standard diff matrix    (np,np)
    real(RDP), intent(in)  :: Jd(:,:,:,:)  !< Jacobian determinant (np,np,np,ne)
    real(RDP), intent(in)  :: G(:,:,:,:,:) !< Laplacian metrics  (np,np,np,ne,6)
    real(RDP), intent(in)  :: lambda       !< Helmholtz parameter
    real(RDP), intent(in)  :: nu           !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:)   !< operand              (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)   !< result               (np,np,np,ne)

    real(RDP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    contiguous :: Ms, Ds, Jd, G, u, v, Ji_n, ub, qb

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_set_1, xmm_set_2, xmm_set_3
    type(LIBXSMM_DMMFunction), save :: xmm_add_1, xmm_add_2, xmm_add_3

    real(RDP), parameter :: ONE  = 1
    real(RDP), parameter :: ZERO = 0

    real(RDP), dimension(size(Ms), size(Ms), size(Ms)) :: M, r, s, t, z

    real(RDP), dimension(size(Ms), size(Ms)) :: Ds_t

    integer   :: np, ne
    integer   :: e, f, i, j, k
    logical   :: get_traces

    !---------------------------------------------------------------------------
    ! initialization

    np = size(Ms)
    ne = size(u,4)

    get_traces = present(Ji_n) .and. present(ub) .and. present(qb)

    ! 3D standard mass matrix
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) =  Ms(k) * Ms(j) * Ms(i)
    end do
    end do
    end do

    ! initialize operators

    Ds_t = transpose(Ds)


    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_set_1, np   , np**2, np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_1, np   , np**2, np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_2, np   , np   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_2, np   , np   , np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_3, np**2, np   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_3, np**2, np   , np, alpha=ONE, beta=ONE )
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      ! v = λ M Jd u ...........................................................

      if (lambda == ZERO) then
        do k = 1, np
        do j = 1, np
        do i = 1, np
          v(i,j,k,e) = ZERO
        end do
        end do
        end do
      else
        do k = 1, np
        do j = 1, np
        do i = 1, np
          v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
        end do
        end do
        end do
      end if

      ! standard derivatives of uᵉ .............................................

      call LIBXSMM_DMMCall(xmm_set_1, Ds, u(:,:,:,e), r(:,:,:))

      do k = 1, np
        call LIBXSMM_DMMCall(xmm_set_2, u(:,:,k,e), Ds_t, s(:,:,k))
      end do

      call LIBXSMM_DMMCall(xmm_set_3, u(:,:,:,e), Ds_t, t(:,:,:))

      ! traces on element boundary .............................................

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

      ! application of metric terms and second differentiation .................
      ! completion: v += ∇ˢ⋅(ν M G ⋅ ∇ˢuᵉ)

      ! direction 1:  z = ν M G(1,:) ⋅ ∇ˢuᵉ ....................................

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
      call LIBXSMM_DMMCall(xmm_add_1, Ds_t, z(:,:,:), v(:,:,:,e))


      ! direction 2: z = ν M G(2,:) ⋅ ∇ˢuᵉ .....................................

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
        call LIBXSMM_DMMCall(xmm_add_2, z(:,:,k), Ds, v(:,:,k,e))
      end do


      ! direction 3: z = ν M G(3,:) ⋅ ∇ˢuᵉ .....................................

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
      call LIBXSMM_DMMCall(xmm_add_3, z(:,:,:), Ds, v(:,:,:,e))


    end do

  end subroutine TPO_Elliptic_DLCI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Elliptic__3D_DLCI__XSMM_RDP
