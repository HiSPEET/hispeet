!> summary:  3d generic element elliptic operator (RLVI)
!> author:   Karl Schoppmann, Joerg Stiller
!> date:     2018/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Elliptic_RLVI_Gen_RWP( Ms, Ds, dx, lambda &
                                    , nu,  u,  v         &
                                    , nub, ub, qb        )

  real(RWP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
  real(RWP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
  real(RWP), intent(in)  :: dx(3)       !< element extensions
  real(RWP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RWP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
  real(RWP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

  real(RWP), optional, intent(inout) :: nub(:,:,:,:) !< ν at element boundary
  real(RWP), optional, intent(inout) :: ub (:,:,:,:) !< u at element boundary
  real(RWP), optional, intent(inout) :: qb (:,:,:,:) !< element boundary flux

  contiguous :: Ms, Ds, nu, u, v, nub, ub, qb

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_nu, z
  real(RWP), dimension(size(Ms), size(Ms))           :: Ds_t

  real(RWP) :: g(3), tmp
  integer   :: np, ne
  integer   :: e, i, j, k, p

  logical :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(u,4)

  get_traces = present(nub) .and. present(ub) .and. present(qb)

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! coefficients
  g = 4 / dx**2

  ! initialize operators
  Ds_t = transpose(Ds)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! vᵉ = λ Mᵉ uᵉ, M_nu = Mᵉνᵉ ................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) = lambda * M(i,j,k) * u(i,j,k,e)
      M_nu(i,j,k) = M(i,j,k) * nu(i,j,k,e)
    end do
    end do
    end do

    ! direction 1 ............................................................

    ! z = [I x I x Dˢ] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds_t(p,i) * u(p,j,k,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    if (get_traces) then
      ! nub = ν;  ub = u;  qb = n⋅ν∇u   @ Γ₁ + Γ₂
      do k = 1, np
      do j = 1, np
        nub(j,k,1,e) =  nu( 1,j,k,e)
        nub(j,k,2,e) =  nu(np,j,k,e)
        ub (j,k,1,e) =  u ( 1,j,k,e)
        ub (j,k,2,e) =  u (np,j,k,e)
        qb (j,k,1,e) = -2/dx(1) * nu( 1,j,k,e) * z( 1,j,k)
        qb (j,k,2,e) =  2/dx(1) * nu(np,j,k,e) * z(np,j,k)
      end do
      end do
    end if

    ! z = g₁ [νᵉMᵉ] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = g(1) * M_nu(i,j,k) * z(i,j,k)
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

    ! direction 2 ............................................................

    ! z = [I x Dˢ x I] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds_t(p,j) * u(i,p,k,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    if (get_traces) then
      ! nub = ν;  ub = u;  qb = n⋅ν∇u   @ Γ₃ + Γ₄
      do k = 1, np
      do i = 1, np
        nub(i,k,3,e) =  nu(i, 1,k,e)
        nub(i,k,4,e) =  nu(i,np,k,e)
        ub (i,k,3,e) =  u (i, 1,k,e)
        ub (i,k,4,e) =  u (i,np,k,e)
        qb (i,k,3,e) = -2/dx(2) * nu(i, 1,k,e) * z(i, 1,k)
        qb (i,k,4,e) =  2/dx(2) * nu(i,np,k,e) * z(i,np,k)
      end do
      end do
    end if

    ! z = g₂ [νᵉMᵉ] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = g(2) * M_nu(i,j,k) * z(i,j,k)
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

    ! direction 3 ............................................................

    ! z = [Dˢ x I x I] uᵉ
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Ds_t(p,k) * u(i,j,p,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    if (get_traces) then
      ! nub = ν;  ub = u;  qb = n⋅ν∇u   @ Γ₅ + Γ₆
      do j = 1, np
      do i = 1, np
        nub(i,j,5,e) =  nu(i,j, 1,e)
        nub(i,j,6,e) =  nu(i,j,np,e)
        ub (i,j,5,e) =  u (i,j, 1,e)
        ub (i,j,6,e) =  u (i,j,np,e)
        qb (i,j,5,e) = -2/dx(3) * nu(i,j, 1,e) * z(i,j, 1)
        qb (i,j,6,e) =  2/dx(3) * nu(i,j,np,e) * z(i,j,np)
      end do
      end do
    end if

    ! z = g₃ [νᵉMᵉ] z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k) = g(3) * M_nu(i,j,k) * z(i,j,k)
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

end subroutine TPO_Elliptic_RLVI_Gen_RWP
