!> summary:   3d generic element diffusion operator (RLCI)
!> author:    Jörg Stiller
!> date:      2019/12/06
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Diffusion_RLCI_Gen_RWP(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
  real(RWP), intent(in)  :: Ms(:)      !< 1D standard mass matrix         (np)
  real(RWP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix (np,np)
  real(RWP), intent(in)  :: dx(3)      !< element extensions
  real(RWP), intent(in)  :: lambda     !< Helmholtz parameter
  real(RWP), intent(in)  :: nu         !< diffusivity
  real(RWP), intent(in)  :: u(:,:,:,:) !< operand                (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:) !< result                 (np,np,np,ne)

  real(RWP), optional, intent(in)    :: Ds(:,:)     !< 1D standard diff matrix
  real(RWP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
  real(RWP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_u
  real(RWP), dimension(size(Ms), size(Ms))           :: Lm
  real(RWP), dimension(size(Ms))                     :: Bs1, Bs2

  real(RWP) :: c(3), tmp, tmp1, tmp2
  integer   :: np, ne
  integer   :: e, i, j, k, p
  logical   :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(u,4)

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, np
  do i = 1, np
    Lm(i,j) = Ls(i,j) / Ms(i)
  end do
  end do

  ! coefficients
  c = 4 * nu / dx**2

  ! handling of boundary traces
  get_traces = present(Ds) .and. present(ub) .and. present(qb)

  if (get_traces) then
    Bs1 = -nu * Ds( 1,:)  ! n⋅νDˢ @ "left" boundary
    Bs2 =  nu * Ds(np,:)  ! n⋅νDˢ @ "right" boundary
  end if

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! M u and lambda M u .......................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M_u(i,j,k) = M(i,j,k) * u(i,j,k,e)
      v(i,j,k,e) = lambda * M_u(i,j,k)
    end do
    end do
    end do

    ! direction 1 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,i) * M_u(p,j,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + c(1) * tmp
    end do
    end do
    end do

    ! direction 2 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,j) * M_u(i,p,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + c(2) * tmp
    end do
    end do
    end do

    ! direction 3 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,k) * M_u(i,j,p)
      end do
      v(i,j,k,e) = v(i,j,k,e) + c(3) * tmp
    end do
    end do
    end do

    ! extract/compute traces ...................................................

    if (get_traces) then

      ! ub = u, n⋅ν∇u  @ Γ₁ + Γ₂
      do k = 1, np
      do j = 1, np
        tmp1 = 0
        tmp2 = 0
        do p = 1, np
          tmp1 = tmp1 + Bs1(p) * u(p,j,k,e)
          tmp2 = tmp2 + Bs2(p) * u(p,j,k,e)
        end do
        ub(j,k,1,e) = u( 1,j,k,e)
        ub(j,k,2,e) = u(np,j,k,e)
        qb(j,k,1,e) = 2/dx(1) * tmp1
        qb(j,k,2,e) = 2/dx(1) * tmp2
      end do
      end do

      ! ub = u, n⋅ν∇u  @ Γ₃ + Γ₄
      do k = 1, np
      do i = 1, np
        tmp1 = 0
        tmp2 = 0
        do p = 1, np
          tmp1 = tmp1 + Bs1(p) * u(i,p,k,e)
          tmp2 = tmp2 + Bs2(p) * u(i,p,k,e)
        end do
        ub(i,k,3,e) = u(i, 1,k,e)
        ub(i,k,4,e) = u(i,np,k,e)
        qb(i,k,3,e) = 2/dx(2) * tmp1
        qb(i,k,4,e) = 2/dx(2) * tmp2
      end do
      end do

      ! ub = u, n⋅ν∇u  @ Γ₅ + Γ₆
      do j = 1, np
      do i = 1, np
        tmp1 = 0
        tmp2 = 0
        do p = 1, np
          tmp1 = tmp1 + Bs1(p) * u(i,j,p,e)
          tmp2 = tmp2 + Bs2(p) * u(i,j,p,e)
        end do
        ub(i,j,5,e) = u(i,j, 1,e)
        ub(i,j,6,e) = u(i,j,np,e)
        qb(i,j,5,e) = 2/dx(3) * tmp1
        qb(i,j,6,e) = 2/dx(3) * tmp2
      end do
      end do

    end if

  end do
  !$omp end do

end subroutine TPO_Diffusion_RLCI_Gen_RWP
