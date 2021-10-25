!-------------------------------------------------------------------------------
!> 3d generic element diffusion operator using hand-crafted suboperators (RLVI)

subroutine PROC(TPO_Diffusion_RLVI_Hand__,_NP_) &
    (ne, Ms, Ds, dx, lambda, nu, u, v, ub, qb)

  integer,   intent(in)  :: ne                    !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)              !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)         !< 1D standard diff matrix
  real(RWP), intent(in)  :: dx(3)                 !< element extensions
  real(RWP), intent(in)  :: lambda                !< Helmholtz parameter
  real(RWP), intent(in)  :: nu(_NP_,_NP_,_NP_,ne) !< diffusivity
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)  !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne)  !< result

  real(RWP), optional, intent(inout) :: ub(_NP_,_NP_,6,ne) !< elem bound values
  real(RWP), optional, intent(inout) :: qb(_NP_,_NP_,6,ne) !< elem bound fluxes

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_nu, z
  real(RWP), dimension(size(Ms), size(Ms))           :: Ds_t

  real(RWP), parameter :: ONE  = 1
  real(RWP), parameter :: ZERO = 0
  real(RWP) :: g(3), tmp

  integer, parameter :: NP_E = _NP_**3
  integer :: e, i, j, k
  logical :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  get_traces = present(ub) .and. present(qb)

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, _NP_
  do j = 1, _NP_
  do i = 1, _NP_
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! coefficients
  g = 4 / dx**2

  ! initialize operators
  Ds_t = transpose(Ds)

  ! make sure that aux array does not contain NaNS
  z = 0

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! vᵉ = λ Mᵉ uᵉ, M_nu = Mᵉνᵉ ................................................

    call SetOperands(lambda, M, nu(:,:,:,e), u(:,:,:,e), v(:,:,:,e), M_nu)

    ! direction 1 ..............................................................

    ! z = [I x I x Dˢ] uᵉ
    call PROC(IxIxQt__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), z)

    if (get_traces) then
      ! ub = u       @ Γ₁ + Γ₂
      ! qb = n⋅ν∇u   @ Γ₁ + Γ₂
      do k = 1, _NP_
      do j = 1, _NP_
        ub(j,k,1,e) =  u(  1 ,j,k,e)
        ub(j,k,2,e) =  u(_NP_,j,k,e)
        qb(j,k,1,e) = -2/dx(1) * nu(  1 ,j,k,e) * z(  1 ,j,k)
        qb(j,k,2,e) =  2/dx(1) * nu(_NP_,j,k,e) * z(_NP_,j,k)
      end do
      end do
    end if

    ! z = g₁ [νᵉMᵉ] z
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = g(1) * M_nu(i,j,k) * z(i,j,k)
    end do
    end do
    end do

    ! vᵉ += [I x I x (Dˢ)ᵀ] z
    call PROC(IxIxQt__,_NP_)(Ds, ONE, ONE, z, v(:,:,:,e))

    ! direction 2 ..............................................................

    ! z = [I x Dˢ x I] uᵉ
    call PROC(IxQtxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), z)

    if (get_traces) then
      ! ub = u       @ Γ₃ + Γ₄
      ! qb = n⋅ν∇u   @ Γ₃ + Γ₄
      do k = 1, _NP_
      do i = 1, _NP_
        ub(i,k,3,e) =  u(i,  1 ,k,e)
        ub(i,k,4,e) =  u(i,_NP_,k,e)
        qb(i,k,3,e) = -2/dx(2) * nu(i,  1 ,k,e) * z(i,  1 ,k)
        qb(i,k,4,e) =  2/dx(2) * nu(i,_NP_,k,e) * z(i,_NP_,k)
      end do
      end do
    end if

    ! z = g₂ [νᵉMᵉ] z
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = g(2) * M_nu(i,j,k) * z(i,j,k)
    end do
    end do
    end do

    ! vᵉ += [I x (Dˢ)ᵀ x I] z
    call PROC(IxQtxI__,_NP_)(Ds, ONE, ONE, z, v(:,:,:,e))

    ! direction 3 ..............................................................

    ! z = [Dˢ x I x I] uᵉ
    call PROC(QtxIxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), z)

    if (get_traces) then
      ! ub = u       @ Γ₅ + Γ₆
      ! qb = n⋅ν∇u   @ Γ₅ + Γ₆
      do j = 1, _NP_
      do i = 1, _NP_
        ub(i,j,5,e) =  u(i,j,  1 ,e)
        ub(i,j,6,e) =  u(i,j,_NP_,e)
        qb(i,j,5,e) = -2/dx(3) * nu(i,j,  1 ,e) * z(i,j,  1 )
        qb(i,j,6,e) =  2/dx(3) * nu(i,j,_NP_,e) * z(i,j,_NP_)
      end do
      end do
    end if

    ! z = g₃ [νᵉMᵉ] z
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = g(3) * M_nu(i,j,k) * z(i,j,k)
    end do
    end do
    end do

    ! vᵉ += [(Dˢ)ᵀ x I x I] z
    call PROC(QtxIxI__,_NP_)(Ds, ONE, ONE, z, v(:,:,:,e))

  end do
  !$omp end do

contains

  !-----------------------------------------------------------------------------
  !> Computation of Helmholtz term and projected diffusivity

  subroutine SetOperands(lambda, M, nu, u, v, M_nu)
    real(RWP), intent(in)  :: lambda
    real(RWP), intent(in)  :: M    (NP_E)
    real(RWP), intent(in)  :: nu   (NP_E)
    real(RWP), intent(in)  :: u    (NP_E)
    real(RWP), intent(out) :: v    (NP_E)
    real(RWP), intent(out) :: M_nu (NP_E)

    integer :: l

    !DIR$ SIMD
    do l = 1, NP_E
      v(l) = lambda * M(l) * u(l)
      M_nu(l) = M(l) * nu(l)
    end do

  end subroutine SetOperands

end subroutine PROC(TPO_Diffusion_RLVI_Hand__,_NP_)
