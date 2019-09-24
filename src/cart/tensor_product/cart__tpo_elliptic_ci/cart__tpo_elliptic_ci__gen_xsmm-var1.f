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

  use Kind_Parameters,  only: RNP
  use ISO_C_Binding,    only: C_Loc  
  use LIBXSMM,          only: LIBXSMM_Dispatch,  LIBXSMM_DMMFunction, &
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
  type(LIBXSMM_DMMFunction) :: xmm_1, xmm_2, xmm_3

  real(RNP), parameter :: ZERO = 0, ONE = 1

  real(RNP) :: M(np,np,np), M_u(np,np,np), Lm(np,np), Lm_t(np,np)
  real(RNP) :: v_1(np,np,np), v_2(np,np,np), v_3(np,np,np)
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
  Lm_t = transpose(Lm)  

  ! coefficients
  g = 4 * nu / dx**2

  ! dispatch LIBXSMM functions 
  call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=ONE, beta=ZERO)
  call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=ONE, beta=ZERO)
  call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=ONE, beta=ZERO)
  
  if ( .not. ( LIBXSMM_Available(xmm_1) .and. &
               LIBXSMM_Available(xmm_2) .and. &
               LIBXSMM_Available(xmm_3) )     ) then
               
		stop "CART__TPO_Elliptic_CI__gen_xsmm: LIBXSMM _Dispatch failed"

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
    end do
    end do
    end do

    ! direction 1 ..............................................................

    !$acc loop collapse(3) independent vector

    call LIBXSMM_MMCall(xmm_1, C_Loc(Lm_t), C_Loc(M_u(1,1,1)), C_Loc(v_1(1,1,1)))

    ! direction 2 ..............................................................
   
    do k = 1, np
      call LIBXSMM_MMCall(xmm_2, C_Loc(M_u(1,1,k)), C_Loc(Lm), C_Loc(v_2(1,1,k))) 
    end do
    
    ! direction 3 ..............................................................
    
    call LIBXSMM_MMCall(xmm_3, C_Loc(M_u(1,1,1)), C_Loc(Lm), C_Loc(v_3(1,1,1))) 
  
    ! assembly of the result ...................................................
      
    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) =                   &
          lambda * M_u(i,j,k)        &      
          + g(1) * v_1(i,j,k)        &
          + g(2) * v_2(i,j,k)        &
          + g(3) * v_3(i,j,k)
    end do
    end do
    end do
      
  end do

!===============================================================================

end subroutine CART__TPO_Elliptic_CI__gen_xsmm
