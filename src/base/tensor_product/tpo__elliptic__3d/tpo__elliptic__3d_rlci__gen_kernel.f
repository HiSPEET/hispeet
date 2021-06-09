!> summary:   3d generic elliptic element operator (RLCI)
!> author:    Jörg Stiller
!> date:      2019/12/06
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Elliptic_RLCI_Gen_RWP(Ms, Ls, lambda, nu, dx, u, v)
  real(RWP), intent(in)  :: Ms(:)      !< 1D standard mass matrix         (np)
  real(RWP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix (np,np)
  real(RWP), intent(in)  :: lambda     !< Helmholtz parameter
  real(RWP), intent(in)  :: nu         !< diffusivity
  real(RWP), intent(in)  :: dx(3)      !< element extensions
  real(RWP), intent(in)  :: u(:,:,:,:) !< operand                (np,np,np,ne)
  real(RWP), intent(out) :: v(:,:,:,:) !< result                 (np,np,np,ne)

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_u
  real(RWP), dimension(size(Ms), size(Ms))           :: Lm

  real(RWP) :: c(3), tmp
  integer   :: np, ne
  integer   :: e, i, j, k, p

  !---------------------------------------------------------------------------
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

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do private(e)
  do e = 1, ne

    ! M u and lambda M u .....................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M_u(i,j,k) = M(i,j,k) * u(i,j,k,e)
      v(i,j,k,e) = lambda * M_u(i,j,k)
    end do
    end do
    end do

    ! direction 1 ............................................................

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

    ! direction 2 ............................................................

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

    ! direction 3 ............................................................

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

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine TPO_Elliptic_RLCI_Gen_RWP
