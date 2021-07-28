!> summary:  3d generic curvilinear element diffusion operator (DLCI)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/07/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Diffusion_DLCI_Gen_RWP(Ms, Ds, Jd, G, lambda, nu, u, v)
  real(RWP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
  real(RWP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
  real(RWP), intent(in)  :: Jd(:,:,:,:) !< Jacobian determinant  (np,np,np,ne)
  real(RWP), intent(in)  :: G(:,:,:,:,:)!< Laplacian metrics     (np,np,np,ne,6)
  real(RWP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RWP), intent(in)  :: nu          !< diffusivity
  real(RWP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, r, s, t, z

  real(RWP) :: tmp
  integer   :: np, ne
  integer   :: e, i, j, k, p

  !---------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(u,4)

  ! element mass matrix

  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! lambda M u .....................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
    end do
    end do
    end do

    ! direction 1: r = ν(MxMxD)u ...............................................

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

    ! direction 2: s = ν(MxDxM)u ...............................................

    !$acc loop collapse(3) independent vector
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

    ! direction 3: t = ν(DxMxM)u ...............................................

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

    ! application of metric terms and second differentiation
    ! direction 1 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,1) * r(i,j,k) &
                                 + G(i,j,k,e,2) * s(i,j,k) &
                                 + G(i,j,k,e,3) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += g₁ [I x I x (MˢDˢ)ᵀ] z
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

    ! direction 2 ............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,2) * r(i,j,k) &
                                 + G(i,j,k,e,4) * s(i,j,k) &
                                 + G(i,j,k,e,5) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += g₂ [I x (MˢDˢ)ᵀ x I] z   ????
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

    ! direction 3 ............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,3) * r(i,j,k) &
                                 + G(i,j,k,e,5) * s(i,j,k) &
                                 + G(i,j,k,e,6) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += g₃ [(MˢDˢ)ᵀ x I x I] z
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
