!-------------------------------------------------------------------------------
!> 3d generic elliptic element operator using hand-crafted suboperators (RLVI)

subroutine PROC(TPO_Elliptic_RLVI_Hand__,_NP_)(ne, Ms, Ds, lambda, nu, dx, u, v)
  integer,   intent(in)  :: ne                    !< num elements
  real(RNP), intent(in)  :: Ms(_NP_)              !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)         !< 1D standard diff matrix
  real(RNP), intent(in)  :: lambda                !< Helmholtz parameter
  real(RNP), intent(in)  :: nu(_NP_,_NP_,_NP_,ne) !< diffusivity
  real(RNP), intent(in)  :: dx(3)                 !< element extensions
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)  !< operand
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne)  !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_u, z
  real(RNP), dimension(size(Ms), size(Ms))           :: Ms_Ds, Msi_Dst

  real(RNP), parameter :: ONE  = 1
  real(RNP), parameter :: ZERO = 0
  real(RNP) :: g(3), tmp

  integer, parameter :: NP_E = _NP_**3
  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, _NP_
  do j = 1, _NP_
  do i = 1, _NP_
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! modified differentiation operators
  do j = 1, _NP_
  do i = 1, _NP_
    Ms_Ds   (i,j) = Ms(i) * Ds(i,j)  ! =  Mˢ     Dˢ    =  [ (Mˢ  Dˢ)ᵀ  ]ᵀ
    Msi_Dst (i,j) = Ds(j,i) / Ms(i)  ! = (Mˢ)⁻¹ (Dˢ)ᵀ  =  [  Dˢ (Mˢ)⁻¹ ]ᵀ
  end do
  end do

  ! coefficients
  g = 4 / dx**2

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(g,M,Ms_Ds,Msi_Dst)
  !$acc parallel
  !$acc loop gang worker private(M_u,z)

  !$omp do private(e)
  do e = 1, ne

    ! M_u = Mᵉuᵉ, vᵉ = λ M uᵉ
    call SetOperands(lambda, M, u(:,:,:,e), M_u, v(:,:,:,e))

    ! vᵉ += g₁ [I x I x (MˢDˢ)ᵀ] [νᵉ] [I x I x (Mˢ)⁻¹(Dˢ)ᵀ)] Mᵉuᵉ
    call PROC(IxIxQt__,_NP_)(Msi_Dst, ONE, ZERO, M_u, z)
    call CompMult(z, nu(:,:,:,e))
    call PROC(IxIxQt__,_NP_)(Ms_Ds, g(1), ONE, z, v(:,:,:,e))

    ! vᵉ += g₂ [I x (MˢDˢ)ᵀ x I] [νᵉ] [I x (Mˢ)⁻¹(Dˢ)ᵀ x I] Mᵉuᵉ
    call PROC(IxQtxI__,_NP_)(Msi_Dst, ONE, ZERO, M_u, z)
    call CompMult(z, nu(:,:,:,e))
    call PROC(IxQtxI__,_NP_)(Ms_Ds, g(2), ONE, z, v(:,:,:,e))

    ! vᵉ += g₃ [(MˢDˢ)ᵀ x I x I] [νᵉ] [(Mˢ)⁻¹(Dˢ)ᵀ x I x I] Mᵉuᵉ
    call PROC(QtxIxI__,_NP_)(Msi_Dst, ONE, ZERO, M_u, z)
    call CompMult(z, nu(:,:,:,e))
    call PROC(QtxIxI__,_NP_)(Ms_Ds, g(3), ONE, z, v(:,:,:,e))

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

contains

  !-----------------------------------------------------------------------------
  !> Projection to element mass matrix and computation of the Helmholtz term

  subroutine SetOperands(lambda, M, u, M_u, v)
    !$acc routine vector
    real(RNP), intent(in)  :: lambda
    real(RNP), intent(in)  :: M   (NP_E)
    real(RNP), intent(in)  :: u   (NP_E)
    real(RNP), intent(out) :: M_u (NP_E)
    real(RNP), intent(out) :: v   (NP_E)

    integer :: l

    !$acc loop vector
    !DIR$ SIMD
    do l = 1, NP_E
      M_u(l) = M(l) * u(l)
      v(l) = lambda * M_u(l)
    end do

  end subroutine SetOperands

  !-----------------------------------------------------------------------------
  !> Component-wise multiplication: z = nu * z

  subroutine CompMult(z, nu)
    !$acc routine vector
    real(RNP), intent(inout) :: z  (NP_E)
    real(RNP), intent(in)    :: nu (NP_E)

    integer :: l

    !$acc loop vector
    !DIR$ SIMD
    do l = 1, NP_E
      z(l) = nu(l) * z(l)
    end do

  end subroutine CompMult

end subroutine PROC(TPO_Elliptic_RLVI_Hand__,_NP_)
