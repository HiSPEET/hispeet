!> summary:  3d generic elliptic element operator (RLVI)
!> author:   Karl Schoppmann, Joerg Stiller
!> date:     2018/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic_3d_RLVI__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Elliptic_RLVI_Gen

contains

  !-----------------------------------------------------------------------------
  !> Generic elliptic element operator (RLVI)

  subroutine TPO_Elliptic_RLVI_Gen(Ms, Ds, lambda, nu, dx, u, v)
    real(RNP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RNP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RNP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
    real(RNP), intent(in)  :: dx(3)       !< element extensions
    real(RNP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RNP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    !---------------------------------------------------------------------------
    ! local variables

    real(RNP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_u, z
    real(RNP), dimension(size(Ms), size(Ms))           :: Ms_Ds, Msi_Dst

    real(RNP) :: g(3), tmp
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

    ! modified differentiation operators
    !! to optimize cache use, the operators are used in transposed form
    !! to save loads, both operators are computed in one loop-nest
    do j = 1, np
    do i = 1, np
      Ms_Ds   (i,j) = Ms(i) * Ds(i,j)  ! = Ms Ds     =  [ (Ms Ds)ᵀ ]ᵀ
      Msi_Dst (i,j) = Ds(j,i) / Ms(i)  ! = Ms⁻¹ Dsᵀ  =  [  Ds Ms⁻¹ ]ᵀ
    end do
    end do

    ! coefficients
    g = 4 / dx**2

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data present(u,v) copyin(g,M,Ms_Ds,Msi_Dst)
    !$acc parallel
    !$acc loop gang worker private(M_u)

    !$omp do private(e)
    do e = 1, ne

      ! M_u = Mᵉuᵉ, vᵉ = λ M uᵉ ................................................

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

      ! z = [νᵉ] [I x I x (Mˢ)⁻¹(Dˢ)ᵀ)] Mᵉuᵉ
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Msi_Dst(p,i) * M_u(p,j,k)
        end do
        z(i,j,k) = nu(i,j,k,e) * tmp
      end do
      end do
      end do

      ! vᵉ += g₁ [I x I x (MˢDˢ)ᵀ] z
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,i) * z(p,j,k)
        end do
        v(i,j,k,e) = v(i,j,k,e) + g(1) * tmp
      end do
      end do
      end do

      ! direction 2 ............................................................

      ! z = [νᵉ] [I x (Mˢ)⁻¹(Dˢ)ᵀ x I] Mᵉuᵉ
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Msi_Dst(p,j) * M_u(i,p,k)
        end do
        z(i,j,k) = nu(i,j,k,e) * tmp
      end do
      end do
      end do

      ! vᵉ += g₂ [I x (MˢDˢ)ᵀ x I] z
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,j) * z(i,p,k)
        end do
        v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp
      end do
      end do
      end do

      ! direction 3 ............................................................

      ! z = [νᵉ] [(Mˢ)⁻¹(Dˢ)ᵀ x I x I] Mᵉuᵉ
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Msi_Dst(p,k) * M_u(i,j,p)
        end do
        z(i,j,k) = nu(i,j,k,e) * tmp
      end do
      end do
      end do

      ! vᵉ += g₃ [(MˢDˢ)ᵀ x I x I] z
      !$acc loop collapse(3) independent vector
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,k) * z(i,j,p)
        end do
        v(i,j,k,e) = v(i,j,k,e) + g(3) * tmp
      end do
      end do
      end do

    end do
    !$omp end do

    !$acc end parallel
    !$acc end data

  end subroutine TPO_Elliptic_RLVI_Gen

  !=============================================================================

end module TPO__Elliptic_3d_RLVI__Gen
