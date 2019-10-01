!> summary:  Elliptic element operator: 3D Cartesian equidistant, LIBXSMM
!> author:   Erik Pfister, Joerg Stiller
!> date:     2019/09/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Elliptic element operator: 3D Cartesian equidistant, LIBXSMM
!>
!> (add documentation here)
!===============================================================================

subroutine CART__TPO_Elliptic_CI__gen_xsmm(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  use LIBXSMM,         only: LIBXSMM_Dispatch,  LIBXSMM_DMMFunction, &
                             LIBXSMM_Available, LIBXSMM_MMCall
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RNP), intent(in)  :: Ms(np)         !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ls(np,np)      !< 1D standard stiffness matrix
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
  real(RNP), intent(in)  :: dx(3)          !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne) !< operand
  real(RNP), intent(out) :: v(np,np,np,ne) !< result

  !-----------------------------------------------------------------------------
  ! local variables

  ! LIBXSMM function pointers
  type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

  real(RNP), parameter :: ZERO = 0, ONE = 1

  real(RNP) :: M(np,np,np), M_u(np,np,np), Lm(np,np)
  real(RNP) :: Lm_1t(np,np), Lm_2(np,np), Lm_3(np,np)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

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
  g = 4 * nu / dx**2

  !Stiffness Matrix directions
  Lm_1t = g(1) * transpose(Lm)
  Lm_2  = g(2) * Lm
  Lm_3  = g(3) * Lm

  ! dispatch LIBXSMM functions 
  !$omp master
  call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=ONE, beta=ONE)
  call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=ONE, beta=ONE)
  call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=ONE, beta=ONE)
  
  if ( .not. ( LIBXSMM_Available(xmm_1) .and. &
               LIBXSMM_Available(xmm_2) .and. &
               LIBXSMM_Available(xmm_3) )     ) then
               
		stop "CART__TPO_Elliptic_CI__gen_xsmm: LIBXSMM _Dispatch failed"

  end if
  !$omp end master
  !$omp barrier
     
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

    call LIBXSMM_DMMCall(xmm_1, Lm_1t, M_u, v(:,:,:,e))

    ! direction 2 ..............................................................
   
    do k = 1, np
      call LIBXSMM_DMMCall(xmm_2, M_u(:,:,k), Lm_2, v(:,:,k,e))) 
    end do
    
    ! direction 3 ..............................................................
    
    call LIBXSMM_DMMCall(xmm_3, M_u, Lm_3, v(:,:,:,e))) 
      
  end do

!===============================================================================

end subroutine CART__TPO_Elliptic_CI__gen_xsmm
