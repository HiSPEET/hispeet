!> summary:  Elliptic element operator: 3D Cartesian equidistant, generic
!> author:   Immo Huismann, Joerg Stiller
!> date:     2015/04/07, revised 2017/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Elliptic element operator: 3D Cartesian equidistant, generic
!>
!> @note
!> When compiled with OpenACC and run on a GPU, any call to this
!> procedure must be followed by an "acc wait" before accessing u
!> or v from CPU or any other than the default accelerator queue.
!> @endnote
!===============================================================================

subroutine CART__TPO_Elliptic_CI__libxsmm_gen(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
  use ISO_C_BINDING,    only: C_LOC  
!#ifdef LIBXSMM 
  use libxsmm,          only: libxsmm_dispatch, libxsmm_dmmfunction, libxsmm_available, libxsmm_mmcall
!#endif
  implicit none

!#ifdef LIBXSMM
 type(libxsmm_dmmfunction) :: xmm_1, xmm_2, xmm_3
!#endif

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

  real(RNP), allocatable :: M(:,:,:), M_u(:,:,:), Lm(:,:)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if(np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! workspace
  allocate(M(np,np,np), M_u(np,np,np), Lm(np,np))

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
    Lm_T=transpose(Lm)

  ! coefficients
  g = 4 * nu / dx**2

  
  call libxsmm_dispatch(xmm_1, np   , np**2, np, &
  &  alpha=ONE, beta=ZERO)
  call libxsmm_dispatch(xmm_2, np   , np   , np, &
  &  alpha=ONE, beta=ZERO)
  call libxsmm_dispatch(xmm_3, np**2, np   , np, &
  &  alpha=ONE, beta=ZERO)
  
  if    (libxsmm_available(xmm_1).and.           &
  &      libxsmm_available(xmm_2).and.           &
  &      libxsmm_available(xmm_3))               &
  &  then
     

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(g,M,Lm) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(M_u)

  !$omp do private(e)
!#ifdef LIBXSMM

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

    call libxsmm_mmcall(xmm_1, &
      & C_LOC(Lm_T), C_LOC(M_u(1,1,1)), C_LOC(v_1(1,1,1)))

    ! direction 2 ..............................................................
   
    do k = 1, np
      call libxsmm_mmcall(xmm_2, &
        & C_LOC(M_u(1,1,k)), C_LOC(Lm), C_LOC(v_2(1,1,k))) 
    end do
    
    ! direction 3 ..............................................................
    
    call libxsmm_mmcall(xmm_3, &
      & C_LOC(M_u(1,1,1)), C_LOC(Lm), C_LOC(v_3(1,1,1))) 
      
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
 
  else
    WRITE(*,*) "Fehler bei libxsmm_dispatch!"
  end if
  
!#endif
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Elliptic_CI__libxsmm_gen
