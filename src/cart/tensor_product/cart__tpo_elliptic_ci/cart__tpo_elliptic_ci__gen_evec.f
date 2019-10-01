!> summary:  Elliptic element operator: 3D Cartesian equidistant, generic EVec
!> author:   Erik Pfister, Joerg Stiller
!> date:     2019/10/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Elliptic element operator: 3D Cartesian equidistant, generic 
!>
!>   * vectorization over elements
!===============================================================================

subroutine CART__TPO_Elliptic_CI__gen_evec(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
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

  integer, parameter :: LVEC    = 4
  integer, parameter :: LVEC_M1 = LVEC - 1
  
  real(RNP) :: M(np,np,np), Lm(np,np)
  real(RNP) :: M_u(0:LVEC_M1,np,np,np)
  real(RNP) :: v_e(0:LVEC_M1,np,np,np)
  real(RNP) :: g(3), tmp
  integer   :: e, i, j, k, l, l_max, p

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

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne, LVEC 
  
    l_max = mod(n1, LVEC)

    ! M u and lambda M u .......................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
       do l = 0, LVEC_M1
          M_u(l,i,j,k) = M(i,j,k) * u(i,j,k,e + min(l,l_max))
          v_e(l,i,j,k) = lambda * M_u(l,i,j,k)
       end do
    end do
    end do
    end do

    ! direction 1 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
!!! ANPASSEN !!!
! Intel SIMD ist klar, aber Vektorisierungsreport erzeugen + anschauen
! GCC vector: lesen, herausfinden wie Vektorisierungsreport erzeugt werden kann
! Alternative für feste np
!   - Template-basierte Parametrisierung, oder
!   - p-Schleife manuell abwickeln (siehe Beispiel von Immo Huismann)
      !DIR$ SIMD
      !GCC$ vector
      do l = 0, LVEC_M1
				tmp = 0
				do p = 1, np
					tmp = tmp + Lm(p,i) * M_u(p,j,k)
				end do
				v_e(l,i,j,k,e) = v_e(l,i,j,k,e) + g(1) * tmp
!!! ANPASSEN ENDE !!!
			end do
    end do
    end do
    end do

    ! direction 2 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
!!! ANPASSEN !!!
      do l = 0, LVEC_M1
				tmp = 0
				do p = 1, np
					tmp = tmp + Lm(p,j) * M_u(i,p,k)
				end do
				v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp
!!! ANPASSEN ENDE !!!
			end do
    end do
    end do
    end do

    ! direction 3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
!!! ANPASSEN !!!
      do l = 0, LVEC_M1
				tmp = 0
				do p = 1, np
					tmp = tmp + Lm(p,k) * M_u(i,j,p)
				end do
				v(i,j,k,e) = v(i,j,k,e) + g(3) * tmp
!!! ANPASSEN ENDE !!!
			end do
    end do
    end do
    end do

    ! assign result ............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      do l = 0, lmax
				v(i,j,k,e+l) = v_e(l,i,j,k)
			end do
    end do
    end do
    end do

  end do

!===============================================================================

end subroutine CART__TPO_Elliptic_CI__gen_evec
