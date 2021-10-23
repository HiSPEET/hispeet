!-------------------------------------------------------------------------------
!> Parametrized 3d diffusion kernel using hand-crafted suboperators (RLCI)

subroutine PROC(TPO_Diffusion_DLCI_Hand__,_NP_) &
    (ne, Ms, Ds, Jd, G, lambda, nu, u, v, Ji_n, ub, qb)

  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)               !< standard mass matrix
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)          !< 1D standard diff matrix
  real(RWP), intent(in)  :: Jd(_NP_,_NP_,_NP_,ne)  !< Jacobian determinant
  real(RWP), intent(in)  :: G(_NP_,_NP_,_NP_,ne,6) !< Laplacian metrics
  real(RWP), intent(in)  :: lambda                 !< Helmholtz parameter λ
  real(RWP), intent(in)  :: nu                     !< diffusivity nu
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne)   !< result

  real(RWP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
  real(RWP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
  real(RWP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

  real(RWP), parameter :: ONE  = 1
  real(RWP), parameter :: ZERO = 0

  real(RWP) :: M(_NP_,_NP_,_NP_), z(_NP_,_NP_,_NP_)
  real(RWP) :: r(_NP_,_NP_,_NP_), s(_NP_,_NP_,_NP_), t(_NP_,_NP_,_NP_)

  real(RWP) :: Ds_t(_NP_,_NP_)

  integer, parameter :: NP_E = _NP_**3
  integer :: e, f, i, j, k
  logical :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  r = 0   ! avoid trouble with NaNs
  s = 0   ! avoid trouble with NaNs
  t = 0   ! avoid trouble with NaNs

  Ds_t = transpose(Ds)

  get_traces = present(Ji_n) .and. present(ub) .and. present(qb)

  ! 3D standard mass matrix
  do k = 1, _NP_
  do j = 1, _NP_
  do i = 1, _NP_
    M(i,j,k) =  Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    ! v = λ M Jd u .............................................................

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
    end do
    end do
    end do

    ! standard derivatives of uᵉ ...............................................

    call PROC(IxIxQt__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), r)
    call PROC(IxQtxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), s)
    call PROC(QtxIxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), t)


    ! traces on element boundary ...............................................

    if (get_traces) then

      ! ub = u,  qb = n⋅ν∇u  @ Γ₁ ∪ Γ₂

      i = 1
      do f = 1, 2
        do k = 1, _NP_
        do j = 1, _NP_
          ub(j,k,f,e) = u(i,j,k,e)
          qb(j,k,f,e) = nu * ( Ji_n(j,k,f,e,1) * r(i,j,k) &
                             + Ji_n(j,k,f,e,2) * s(i,j,k) &
                             + Ji_n(j,k,f,e,3) * t(i,j,k) )
        end do
        end do
        i = _NP_
      end do

      ! ub = u,  qb = n⋅ν∇u  @ Γ₃ ∪ Γ₄

      j = 1
      do f = 3, 4
        do k = 1, _NP_
        do i = 1, _NP_
          ub(i,k,f,e) = u(i,j,k,e)
          qb(i,k,f,e) = nu * ( Ji_n(i,k,f,e,1) * r(i,j,k) &
                             + Ji_n(i,k,f,e,2) * s(i,j,k) &
                             + Ji_n(i,k,f,e,3) * t(i,j,k) )
        end do
        end do
        j = _NP_
      end do

      ! ub = u,  qb = n⋅ν∇u  @ Γ₅ ∪ Γ₆

      k = 1
      do f = 5, 6
        do j = 1, _NP_
        do i = 1, _NP_
          ub(i,j,f,e) = u(i,j,k,e)
          qb(i,j,f,e) = nu * ( Ji_n(i,j,f,e,1) * r(i,j,k) &
                             + Ji_n(i,j,f,e,2) * s(i,j,k) &
                             + Ji_n(i,j,f,e,3) * t(i,j,k) )
        end do
        end do
        k = _NP_
      end do

    end if

    ! application of metric terms and second differentiation ...................
    ! completion: v += ∇ˢ⋅(ν M G ⋅ ∇ˢuᵉ)

    ! direction 1:  z = ν M G(1,:) ⋅ ∇ˢuᵉ ......................................

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,1) * r(i,j,k) &
                                 + G(i,j,k,e,2) * s(i,j,k) &
                                 + G(i,j,k,e,3) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [I x I x (Dˢ)ᵀ] z
    call PROC(IxIxQt__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))


    ! direction 2: z = ν M G(2,:) ⋅ ∇ˢuᵉ .......................................

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,2) * r(i,j,k) &
                                 + G(i,j,k,e,4) * s(i,j,k) &
                                 + G(i,j,k,e,5) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [I x (Dˢ)ᵀ x I] z
    call PROC(IxQtxI__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))


    ! direction 3: z = ν M G(3,:) ⋅ ∇ˢuᵉ .......................................

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,3) * r(i,j,k) &
                                 + G(i,j,k,e,5) * s(i,j,k) &
                                 + G(i,j,k,e,6) * t(i,j,k) )
    end do
    end do
    end do

    ! vᵉ += [(Dˢ)ᵀ x I x I] z
    call PROC(QtxIxI__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))

  end do

end subroutine PROC(TPO_Diffusion_DLCI_Hand__,_NP_)
