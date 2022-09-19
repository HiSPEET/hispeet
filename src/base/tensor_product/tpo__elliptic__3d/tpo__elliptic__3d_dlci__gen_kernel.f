!> summary:  3d generic curvilinear element elliptic operator (DLCI)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/07/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Elliptic_DLCI_Gen_RWP(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                     Ji_n, ub, qb)

  real(RWP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
  real(RWP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
  real(RWP), intent(in)  :: Jd(:,:,:,:) !< Jacobian determinant  (np,np,np,ne)
  real(RWP), intent(in)  :: G(:,:,:,:,:)!< Laplacian metrics   (np,np,np,ne,6)
  real(RWP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RWP), intent(in)  :: nu          !< diffusivity
  real(RWP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

  real(RWP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
  real(RWP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
  real(RWP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

  contiguous :: Ms, Ds, Jd, G, u, v, Ji_n, ub, qb

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: ZERO = 0

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, r, s, t, z

  real(RWP) :: tmp
  integer   :: np, ne
  integer   :: e, f, i, j, k, p
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
    M(i,j,k) = Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! v = λ M Jd u .............................................................

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

    ! standard derivatives of  uᵉ ..............................................

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

    ! application of metric terms and second differentiation ...................
    ! completion: v += ∇ˢ⋅(ν M G ⋅ ∇ˢuᵉ)

    ! direction 1:  z = ν M G(1,:) ⋅ ∇ˢuᵉ ......................................

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

    ! direction 2: z = ν M G(2,:) ⋅ ∇ˢuᵉ .......................................

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

    ! direction 3: z = ν M G(3,:) ⋅ ∇ˢuᵉ .......................................

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

end subroutine TPO_Elliptic_DLCI_Gen_RWP
