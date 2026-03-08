!> summary:   3D element convection operator using LIBXSMM for real(RDP) (D)
!> author:    Jerome Michel, Jörg Stiller
!> date:      2021/01/21 - 2026/03/08
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__INS_Convection__3D_D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_INS_Convection_D_XSMM

  interface TPO_INS_Convection_D_XSMM
    module procedure TPO_INS_Convection_D_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D convection element operator for INS based on LIBXSMM (D)

  subroutine TPO_INS_Convection_D_XSMM_RDP( nv, nq, ne, bc, alpha, D_v, I_vq  &
                                          , w_q , Jd_q, Ji_q, a_q, n_q, v, vp &
                                          , F_c                               )

    integer, intent(in) :: nv !< number of velocity points per element direction
    integer, intent(in) :: nq !< number of quadrature points per direction
    integer, intent(in) :: ne !< number of elements

    character, intent(in) :: bc(6,ne) !< element face BC

    real(RDP), intent(in) :: alpha !< 0/0.5/1: flux/skew-sym/convective form

    real(RDP), intent(in) :: D_v (nv,nv) !< 1D diff operator on v-points
    real(RDP), intent(in) :: I_vq(nq,nv) !< 1D interpolation from to q-points
    real(RDP), intent(in) :: w_q (nq)    !< 1D quadrature weights

    real(RDP), intent(in)  :: Jd_q(nq,nq,nq,ne)     !< Jacobian determinant
    real(RDP), intent(in)  :: Ji_q(nq,nq,nq,ne,3,3) !< inverse Jacobian matrix
    real(RDP), intent(in)  :: a_q(nq,nq,6,ne)       !< face area coefficients
    real(RDP), intent(in)  :: n_q(nq,nq,6,ne,3)     !< face unit normal vectors
    real(RDP), intent(in)  :: v(nv,nv,nv,ne,3)      !< velocity
    real(RDP), intent(in)  :: vp(nv,nv,6,ne,3)      !< velocity v⁺ at faces
    real(RDP), intent(out) :: F_c(nv,nv,nv,ne,3)    !< convective flux integrals

    ! local variables ............................................................

    logical   :: flux_form        ! T for flux form (α=0), F else
    logical   :: interpolate      ! T/F if interpolation is/not required
    logical   :: lobatto          ! T/F if quadrature is/not of Lobatto type

    real(RDP) :: D_vq(nq,nv)      ! interpolated differentiation operator
    real(RDP) :: D_vq_t(nv,nq)    ! transposed interp. different. operator
    real(RDP) :: I_vq_t(nv,nq)    ! transposed interpolation operator
    real(RDP) :: M(nq,nq,nq)      ! mass matrix based on quadrature points

    real(RDP) :: v_qe(nq,nq,nq,3) ! velocity at element quadrature points
    real(RDP) :: vm_qf(nq,nq,3)   ! interior velocity v⁻ at face quadrature points
    real(RDP) :: vp_qf(nq,nq,3)   ! exterior velocity v⁺ at face quadrature points

    real(RDP) :: M_div_v(nq,nq,nq)     ! weighted velocity divergence,
    real(RDP) :: grad_v(nq,nq,nq,3,3)  ! standard velocity gradient, and
    real(RDP) :: M_Ji_vv(nq,nq,nq,3,3) ! MJ⁻¹⋅vv at element quadrature points

    real(RDP) :: y1(nq,nv,nv,2)
    real(RDP) :: y2(nq,nq,nv,3)
    real(RDP) :: z1(nq,nv,nv)
    real(RDP) :: z2(nq,nq,nv)
    real(RDP) :: ff(nq,nq,3)
    real(RDP) :: zd(nv,nv,nv)

    real(RDP) :: v1v1, v1v2, v1v3, v2v2, v2v3, v3v3
    real(RDP) :: vm_1, vm_2, vm_3, vm_n
    real(RDP) :: vp_1, vp_2, vp_3, vp_n
    real(RDP) :: vb_1, vb_2, vb_3, vb_n
    real(RDP) :: cf, vn

    integer   :: c, e, f, i, j, k
    integer   :: iv, iq, jv, jq, kv, kq

    type(LIBXSMM_DMMFunction), save :: xmm_set_1, xmm_set_2, xmm_set_3
    type(LIBXSMM_DMMFunction), save :: xmm_add_1, xmm_add_2, xmm_add_3
    type(LIBXSMM_DMMFunction), save :: xmm_set_1t, xmm_set_2t, xmm_set_3t
    type(LIBXSMM_DMMFunction), save :: xmm_set_f1, xmm_set_f2, xmm_set_f3
    type(LIBXSMM_DMMFunction), save :: xmm_set_f
    type(LIBXSMM_DMMFunction), save :: xmm_add_f

    real(RDP), parameter :: ONE  = 1
    real(RDP), parameter :: ZERO = 0

    !-----------------------------------------------------------------------------
    ! Initialization

    flux_form = alpha == 0
    lobatto = I_vq(1,1) == 1 .and. I_vq(nq,nv) == 1
    interpolate = nq /= nv  .or. .not. lobatto

    if (interpolate) then
      D_vq   = matmul(I_vq, D_v)
      D_vq_t = transpose(D_vq)
      I_vq_t = transpose(I_vq)
    else
      D_vq   = D_v
      D_vq_t = transpose(D_v)
    end if

    if (flux_form) then
      M_div_v = 0
    end if

    ! XSMM dispatch
    !$omp master

    ! for volume integral
    call LIBXSMM_Dispatch(xmm_set_1, nv   , nv**2, nq, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_2, nq   , nv   , nq, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_3, nq**2, nv   , nq, alpha=ONE, beta=ZERO)

    call LIBXSMM_Dispatch(xmm_add_1, nv   , nv**2, nq, alpha=ONE, beta=ONE)
    call LIBXSMM_Dispatch(xmm_add_2, nq   , nv   , nq, alpha=ONE, beta=ONE)
    call LIBXSMM_Dispatch(xmm_add_3, nq**2, nv   , nq, alpha=ONE, beta=ONE)

    call LIBXSMM_Dispatch(xmm_set_1t, nq   , nv**2, nv, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_2t, nq   , nq   , nv, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_3t, nq**2, nq   , nv, alpha=ONE, beta=ZERO)

    ! for face integral
    call LIBXSMM_Dispatch(xmm_set_f1, nq, nv, nv, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_f2, nq, nq, nv, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_f3, nq, nv, nq, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_set_f , nv, nv, nq, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_f , nv, nv, nq, alpha=ONE, beta=ONE)

    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------

    !$omp do
    Elements: do e = 1, ne

      !-------------------------------------------------------------------------
      ! volume integral

      ! v at element quadrature points :::::::::::::::::::::::::::::::::::::::::

      if (interpolate) then

        do c = 1, 3

          ! direction 1
          call LIBXSMM_DMMCall(xmm_set_1t, I_vq, v(:,:,:,e,c), z1)

          ! direction 2
          do k = 1, nv
            call LIBXSMM_DMMCall(xmm_set_2t, z1(:,:,k), I_vq_t, z2(:,:,k))
          end do

          ! direction 3
          call LIBXSMM_DMMCall(xmm_set_3t, z2, I_vq_t, v_qe(:,:,:,c))

        end do

      else

        do c = 1, 3
          do k = 1, nq
          do j = 1, nq
          do i = 1, nq
            v_qe(i,j,k,c) = v(i,j,k,e,c)
          end do
          end do
          end do
        end do

      end if

      ! mass matrix ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      do k = 1, nq
      do j = 1, nq
      do i = 1, nq
        M(i,j,k) = w_q(k) * w_q(j) * w_q(i) * Jd_q(i,j,k,e)
      end do
      end do
      end do

      ! flux contribution ::::::::::::::::::::::::::::::::::::::::::::::::::::::

      ! MJ⁻¹⋅vv ................................................................

      do k = 1, nq
      do j = 1, nq
      do i = 1, nq

        v1v1 = v_qe(i,j,k,1) * v_qe(i,j,k,1)
        v1v2 = v_qe(i,j,k,1) * v_qe(i,j,k,2)
        v1v3 = v_qe(i,j,k,1) * v_qe(i,j,k,3)
        v2v2 = v_qe(i,j,k,2) * v_qe(i,j,k,2)
        v2v3 = v_qe(i,j,k,2) * v_qe(i,j,k,3)
        v3v3 = v_qe(i,j,k,3) * v_qe(i,j,k,3)

        M_Ji_vv(i,j,k,1,1) = M(i,j,k) * ( Ji_q(i,j,k,e,1,1) * v1v1 &
                                        + Ji_q(i,j,k,e,1,2) * v1v2 & ! = v2v1
                                        + Ji_q(i,j,k,e,1,3) * v1v3 ) ! = v3v1

        M_Ji_vv(i,j,k,1,2) = M(i,j,k) * ( Ji_q(i,j,k,e,1,1) * v1v2 &
                                        + Ji_q(i,j,k,e,1,2) * v2v2 &
                                        + Ji_q(i,j,k,e,1,3) * v2v3 ) ! = v3v2

        M_Ji_vv(i,j,k,1,3) = M(i,j,k) * ( Ji_q(i,j,k,e,1,1) * v1v3 &
                                        + Ji_q(i,j,k,e,1,2) * v2v3 &
                                        + Ji_q(i,j,k,e,1,3) * v3v3 )

        M_Ji_vv(i,j,k,2,1) = M(i,j,k) * ( Ji_q(i,j,k,e,2,1) * v1v1 &
                                        + Ji_q(i,j,k,e,2,2) * v1v2 & ! = v2v1
                                        + Ji_q(i,j,k,e,2,3) * v1v3 ) ! = v3v1

        M_Ji_vv(i,j,k,2,2) = M(i,j,k) * ( Ji_q(i,j,k,e,2,1) * v1v2 &
                                        + Ji_q(i,j,k,e,2,2) * v2v2 &
                                        + Ji_q(i,j,k,e,2,3) * v2v3 ) ! = v3v2

        M_Ji_vv(i,j,k,2,3) = M(i,j,k) * ( Ji_q(i,j,k,e,2,1) * v1v3 &
                                        + Ji_q(i,j,k,e,2,2) * v2v3 &
                                        + Ji_q(i,j,k,e,2,3) * v3v3 )

        M_Ji_vv(i,j,k,3,1) = M(i,j,k) * ( Ji_q(i,j,k,e,3,1) * v1v1 &
                                        + Ji_q(i,j,k,e,3,2) * v1v2 & ! = v2v1
                                        + Ji_q(i,j,k,e,3,3) * v1v3 ) ! = v3v1

        M_Ji_vv(i,j,k,3,2) = M(i,j,k) * ( Ji_q(i,j,k,e,3,1) * v1v2 &
                                        + Ji_q(i,j,k,e,3,2) * v2v2 &
                                        + Ji_q(i,j,k,e,3,3) * v2v3 ) ! = v3v2

        M_Ji_vv(i,j,k,3,3) = M(i,j,k) * ( Ji_q(i,j,k,e,3,1) * v1v3 &
                                        + Ji_q(i,j,k,e,3,2) * v2v3 &
                                        + Ji_q(i,j,k,e,3,3) * v3v3 )

      end do
      end do
      end do

      if (interpolate) then

        ! F_c = (I₃×I₂×D₁)ᵀ MJ⁻¹⋅vv ............................................

        do c = 1, 3

          ! note: pass through directions in reverse order to reuse workspace
          ! projection in direction 3
          call LIBXSMM_DMMCall(xmm_set_3, M_Ji_vv(:,:,:,1,c), I_vq, z2)

          ! projection in direction 2
          do k = 1, nv
            call LIBXSMM_DMMCall(xmm_set_2, z2(:,:,k), I_vq, z1(:,:,k))
          end do

          ! differentiation in direction 1
          call LIBXSMM_DMMCall(xmm_set_1, D_vq_t, z1, F_c(:,:,:,e,c))

        end do

        ! F_c += (I₃×D₂×I₁)ᵀ MJ⁻¹⋅vv ...........................................

        do c = 1, 3

          ! projection in direction 3
          call LIBXSMM_DMMCall(xmm_set_3, M_Ji_vv(:,:,:,2,c), I_vq, z2)

          ! differentiation in direction 2
          do k = 1, nv
            call LIBXSMM_DMMCall(xmm_set_2, z2(:,:,k), D_vq, z1(:,:,k))
          end do

          ! projection in direction 1
          call LIBXSMM_DMMCall(xmm_add_1, I_vq_t, z1, F_c(:,:,:,e,c))

        end do

        ! F_c += (D₃×I₂×I₁)ᵀ MJ⁻¹⋅vv ...........................................

        do c = 1, 3

          ! differentiation in direction 3
          call LIBXSMM_DMMCall(xmm_set_3, M_Ji_vv(:,:,:,3,c), D_vq, z2)

          ! projection in direction 2
          do k = 1, nv
            call LIBXSMM_DMMCall(xmm_set_2, z2(:,:,k), I_vq, z1(:,:,k))
          end do

          ! projection in direction 1
          call LIBXSMM_DMMCall(xmm_add_1, I_vq_t, z1, F_c(:,:,:,e,c))

        end do

      else

        ! standard divergence of MJ⁻¹⋅vv (on-grid) .............................

        do  c = 1, 3

          ! F_c = (I₃×I₂×D₁)ᵀ MJ⁻¹⋅vv
          call LIBXSMM_DMMCall(xmm_set_1, D_vq_t, M_Ji_vv(:,:,:,1,c), zd)

          ! F_c += (I₃×D₂×I₁)ᵀ MJ⁻¹⋅vv
          do k = 1, nv
            call LIBXSMM_DMMCall(xmm_add_2, M_Ji_vv(:,:,k,2,c), D_vq, zd(:,:,k))
          end do

          ! F_c += (D₃×I₂×I₁)ᵀ MJ⁻¹⋅vv
          call LIBXSMM_DMMCall(xmm_add_3, M_Ji_vv(:,:,:,3,c), D_vq, zd)

          F_c(:,:,:,e,c) = zd

        end do

      end if

      ! divergence contribution ::::::::::::::::::::::::::::::::::::::::::::::::

      DIVERGENCE: if (.not. flux_form) then

        ! standard gradient of velocity ........................................

        if (interpolate) then

          do c = 1, 3

            ! direction 1
            call LIBXSMM_DMMCall(xmm_set_1t, D_vq, v(:,:,:,e,c), y1(:,:,:,1))
            call LIBXSMM_DMMCall(xmm_set_1t, I_vq, v(:,:,:,e,c), y1(:,:,:,2))

            ! direction 2
            do k = 1, nv
              call LIBXSMM_DMMCall(xmm_set_2t, y1(:,:,k,1), I_vq_t, y2(:,:,k,1))
              call LIBXSMM_DMMCall(xmm_set_2t, y1(:,:,k,2), D_vq_t, y2(:,:,k,2))
              call LIBXSMM_DMMCall(xmm_set_2t, y1(:,:,k,2), I_vq_t, y2(:,:,k,3))
            end do

            ! direction 3
            call LIBXSMM_DMMCall(xmm_set_3t, y2(:,:,:,1), I_vq_t, grad_v(:,:,:,1,c))
            call LIBXSMM_DMMCall(xmm_set_3t, y2(:,:,:,2), I_vq_t, grad_v(:,:,:,2,c))
            call LIBXSMM_DMMCall(xmm_set_3t, y2(:,:,:,3), D_vq_t, grad_v(:,:,:,3,c))

          end do

        else

          do c = 1, 3

            ! direction 1
            call LIBXSMM_DMMCall(xmm_set_1, D_vq, v(:,:,:,e,c), grad_v(:,:,:,1,c))

            ! direction 2
            do k = 1, nv
              call LIBXSMM_DMMCall(xmm_set_2, v(:,:,k,e,c), D_vq_t, grad_v(:,:,k,2,c))
            end do

            ! direction 3
            call LIBXSMM_DMMCall(xmm_set_3, v(:,:,:,e,c), D_vq_t, grad_v(:,:,:,3,c))

          end do

        end if

        ! weighted divergence ..................................................

        do k = 1, nq
        do j = 1, nq
        do i = 1, nq
          M_div_v(i,j,k) =  M(i,j,k) * ( Ji_q(i,j,k,e,1,1) * grad_v(i,j,k,1,1) &
                                       + Ji_q(i,j,k,e,1,2) * grad_v(i,j,k,2,1) &
                                       + Ji_q(i,j,k,e,1,3) * grad_v(i,j,k,3,1) &
                                       + Ji_q(i,j,k,e,2,1) * grad_v(i,j,k,1,2) &
                                       + Ji_q(i,j,k,e,2,2) * grad_v(i,j,k,2,2) &
                                       + Ji_q(i,j,k,e,2,3) * grad_v(i,j,k,3,2) &
                                       + Ji_q(i,j,k,e,3,1) * grad_v(i,j,k,1,3) &
                                       + Ji_q(i,j,k,e,3,2) * grad_v(i,j,k,2,3) &
                                       + Ji_q(i,j,k,e,3,3) * grad_v(i,j,k,3,3) )
        end do
        end do
        end do

        ! F_c += (I₃×I₂×I₁)ᵀ M α(∇⋅v) v_c ......................................

        if (interpolate) then

          associate(w => grad_v(:,:,:,3,1))
            do c = 1, 3

              ! w = M α(∇⋅v) v_c
              w = alpha * M_div_v * v_qe(:,:,:,c)

              ! projection in direction 3
              call LIBXSMM_DMMCall(xmm_set_3, w, I_vq, z2)

              ! projection in direction 2
              do k = 1, nv
                call LIBXSMM_DMMCall(xmm_set_2, z2(:,:,k), I_vq, z1(:,:,k))
              end do

              ! projection in direction 1
              call LIBXSMM_DMMCall(xmm_add_1, I_vq_t, z1, F_c(:,:,:,e,c))

            end do
          end associate

        else

          do k = 1, nv
          do j = 1, nv
          do i = 1, nv
            F_c(i,j,k,e,1) = F_c(i,j,k,e,1) + alpha * M_div_v(i,j,k) * v(i,j,k,e,1)
            F_c(i,j,k,e,2) = F_c(i,j,k,e,2) + alpha * M_div_v(i,j,k) * v(i,j,k,e,2)
            F_c(i,j,k,e,3) = F_c(i,j,k,e,3) + alpha * M_div_v(i,j,k) * v(i,j,k,e,3)
          end do
          end do
          end do

        end if

      end if DIVERGENCE

      !-------------------------------------------------------------------------
      ! face integrals

      ! Γ₁ ∪ Γ₂ ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Faces_1_and_2: do f = 1, 2

        iv = 1 + (f - 1) * (nv - 1)   ! i-index of the face in v-space
        iq = 1 + (f - 1) * (nq - 1)   ! i-index of the face in q-space

        if (interpolate) then
          associate(zf => z1(:,:,1))

            ! v⁻ at q-points ...................................................

            if (lobatto) then

              ! extraction from element q-points
              do c = 1, 3
                do k = 1, nq
                do j = 1, nq
                  vm_qf(j,k,c) = v_qe(iq,j,k,c)
                end do
                end do
              end do

            else

              ! interpolation from v-points to face q-points
              do c = 1, 3

                ! direction 2
                call LIBXSMM_DMMCall(xmm_set_f1, I_vq, v(iv,:,:,e,c), zf(:,:))
                ! direction 3
                call LIBXSMM_DMMCall(xmm_set_f2, zf(:,:), I_vq_t, vm_qf(:,:,c))

              end do

            end if

            ! v⁺ at q-points ...................................................

            ! interpolation from v-points to face q-points
            do c = 1, 3

              ! direction 2
              call LIBXSMM_DMMCall(xmm_set_f1, I_vq, vp(:,:,f,e,c), zf(:,:))
              ! direction 3
              call LIBXSMM_DMMCall(xmm_set_f2, zf(:,:), I_vq_t, vp_qf(:,:,c))

            end do

          end associate

        else ! no interpolation

          ! extract v⁻ and v⁺ ..................................................

          do c = 1, 3
            do k = 1, nv
            do j = 1, nv
              vm_qf(j,k,c) = v_qe(iv,j,k,c)
            end do
            end do
          end do

          do c = 1, 3
            do k = 1, nv
            do j = 1, nv
              vp_qf(j,k,c) = vp(j,k,f,e,c)
            end do
            end do
          end do

        end if

        ! integrand ............................................................

        select case(bc(f,e))

        case('D')

          ! boundary flux, assuming vb = {v} . . . . . . . . . . . . . . . . . . .

          do k = 1, nq
          do j = 1, nq

            vb_1 = 0.5 * (vm_qf(j,k,1) + vp_qf(j,k,1))
            vb_2 = 0.5 * (vm_qf(j,k,2) + vp_qf(j,k,2))
            vb_3 = 0.5 * (vm_qf(j,k,3) + vp_qf(j,k,3))

            vb_n = n_q(j,k,f,e,1) * vb_1  &
                 + n_q(j,k,f,e,2) * vb_2  &
                 + n_q(j,k,f,e,3) * vb_3

            cf = -w_q(j) * w_q(k) * a_q(j,k,f,e) * vb_n

            ff(j,k,1) = cf * vb_1
            ff(j,k,2) = cf * vb_2
            ff(j,k,3) = cf * vb_3

          end do
          end do

        case default

          ! local Lax-Friedrichs flux  . . . . . . . . . . . . . . . . . . . . . .

          do k = 1, nq
          do j = 1, nq

            vm_1 = vm_qf(j,k,1)
            vm_2 = vm_qf(j,k,2)
            vm_3 = vm_qf(j,k,3)

            vm_n = n_q(j,k,f,e,1) * vm_1  &
                 + n_q(j,k,f,e,2) * vm_2  &
                 + n_q(j,k,f,e,3) * vm_3

            vp_1 = vp_qf(j,k,1)
            vp_2 = vp_qf(j,k,2)
            vp_3 = vp_qf(j,k,3)

            vp_n = n_q(j,k,f,e,1) * vp_1  &
                 + n_q(j,k,f,e,2) * vp_2  &
                 + n_q(j,k,f,e,3) * vp_3

            cf = -w_q(j) * w_q(k) * a_q(j,k,f,e)
            vn = max(abs(vm_n), abs(vp_n))

            ff(j,k,1) = cf * (0.5 * (vm_n * vm_1 + vp_n * vp_1) + vn * (vm_1 - vp_1))
            ff(j,k,2) = cf * (0.5 * (vm_n * vm_2 + vp_n * vp_2) + vn * (vm_2 - vp_2))
            ff(j,k,3) = cf * (0.5 * (vm_n * vm_3 + vp_n * vp_3) + vn * (vm_3 - vp_3))

          end do
          end do

        end select

        ! face integral contributions ..........................................

        if (interpolate) then

          associate(zf => z1(:,:,1), F_cf => zd(:,:,1))
            do c = 1, 3

              ! integration over direction 3
              call LIBXSMM_DMMCall(xmm_set_f3, ff(:,:,c), I_vq, zf)

              ! integration over direction 2
              call LIBXSMM_DMMCall(xmm_set_f, I_vq_t, zf, F_cf)
              F_c(iv,:,:,e,c) = F_c(iv,:,:,e,c) + F_cf

            end do
          end associate

        else

          do c = 1, 3
            do k = 1, nv
            do j = 1, nv
              F_c(iv,j,k,e,c) = F_c(iv,j,k,e,c) + ff(j,k,c)
            end do
            end do
          end do

        end if

      end do Faces_1_and_2

      ! Γ₃ ∪ Γ₄ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Faces_3_and_4: do f = 3, 4

        jv = 1 + (f - 3) * (nv - 1)   ! j-index of the face in v-space
        jq = 1 + (f - 3) * (nq - 1)   ! j-index of the face in q-space

        if (interpolate) then
          associate(zf => z1(:,:,1))

            ! v⁻ at q-points ...................................................

            if (lobatto) then

              ! extraction from element q-points
              do c = 1, 3
                do k = 1, nq
                do i = 1, nq
                  vm_qf(i,k,c) = v_qe(i,jq,k,c)
                end do
                end do
              end do

            else

              ! interpolation from v-points to face q-points
              do c = 1, 3

                ! direction 1
                call LIBXSMM_DMMCall(xmm_set_f1, I_vq, v(:,jv,:,e,c), zf)

                ! direction 3
                call LIBXSMM_DMMCall(xmm_set_f2, zf, I_vq_t, vm_qf(:,:,c))

              end do

            end if

            ! v⁺ at q-points ...................................................

            ! interpolation from v-points to face q-points
            do c = 1, 3

              ! direction 1
              call LIBXSMM_DMMCall(xmm_set_f1, I_vq, vp(:,:,f,e,c), zf)

              ! direction 3
              call LIBXSMM_DMMCall(xmm_set_f2, zf, I_vq_t, vp_qf(:,:,c))

            end do

          end associate

        else ! no interpolation

          ! extract v⁻ and v⁺ ..................................................

          do c = 1, 3
            do k = 1, nv
            do i = 1, nv
              vm_qf(i,k,c) = v_qe(i,jv,k,c)
            end do
            end do
          end do

          do c = 1, 3
            do k = 1, nv
            do i = 1, nv
              vp_qf(i,k,c) = vp(i,k,f,e,c)
            end do
            end do
          end do

        end if

        ! integrand ............................................................

        select case(bc(f,e))

        case('D')

          ! boundary flux, assuming vb = {v} . . . . . . . . . . . . . . . . . .

          do k = 1, nq
          do i = 1, nq

            vb_1 = 0.5 * (vm_qf(i,k,1) + vp_qf(i,k,1))
            vb_2 = 0.5 * (vm_qf(i,k,2) + vp_qf(i,k,2))
            vb_3 = 0.5 * (vm_qf(i,k,3) + vp_qf(i,k,3))

            vb_n = n_q(i,k,f,e,1) * vb_1  &
                 + n_q(i,k,f,e,2) * vb_2  &
                 + n_q(i,k,f,e,3) * vb_3

            cf = -w_q(i) * w_q(k) * a_q(i,k,f,e) * vb_n

            ff(i,k,1) = cf * vb_1
            ff(i,k,2) = cf * vb_2
            ff(i,k,3) = cf * vb_3

          end do
          end do

        case default

          ! local Lax-Friedrichs flux  . . . . . . . . . . . . . . . . . . . . .

          do k = 1, nq
          do i = 1, nq

            vm_1 = vm_qf(i,k,1)
            vm_2 = vm_qf(i,k,2)
            vm_3 = vm_qf(i,k,3)

            vm_n = n_q(i,k,f,e,1) * vm_1  &
                 + n_q(i,k,f,e,2) * vm_2  &
                 + n_q(i,k,f,e,3) * vm_3

            vp_1 = vp_qf(i,k,1)
            vp_2 = vp_qf(i,k,2)
            vp_3 = vp_qf(i,k,3)

            vp_n = n_q(i,k,f,e,1) * vp_1  &
                 + n_q(i,k,f,e,2) * vp_2  &
                 + n_q(i,k,f,e,3) * vp_3

            cf = -w_q(i) * w_q(k) * a_q(i,k,f,e)
            vn = max(abs(vm_n), abs(vp_n))

            ff(i,k,1) = cf * (0.5 * (vm_n * vm_1 + vp_n * vp_1) + vn * (vm_1 - vp_1))
            ff(i,k,2) = cf * (0.5 * (vm_n * vm_2 + vp_n * vp_2) + vn * (vm_2 - vp_2))
            ff(i,k,3) = cf * (0.5 * (vm_n * vm_3 + vp_n * vp_3) + vn * (vm_3 - vp_3))

          end do
          end do

        end select

        ! face integral contributions ..........................................

        if (interpolate) then

          associate(zf => z1(:,:,1), F_cf => zd(:,:,1))
            do c = 1, 3

              ! integration over direction 3
              call LIBXSMM_DMMCall(xmm_set_f3, ff(:,:,c), I_vq, zf)

              ! integration over direction 1
              call LIBXSMM_DMMCall(xmm_set_f, I_vq_t, zf, F_cf)
              F_c(:,jv,:,e,c) = F_c(:,jv,:,e,c) + F_cf

            end do
          end associate

        else

          do c = 1, 3
            do k = 1, nv
            do i = 1, nv
              F_c(i,jv,k,e,c) = F_c(i,jv,k,e,c) + ff(i,k,c)
            end do
            end do
          end do

        end if

      end do Faces_3_and_4

      ! Γ₅ ∪ Γ₆ ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Faces_5_and_6: do f = 5, 6

        kv = 1 + (f - 5) * (nv - 1)   ! k-index of the face in v-space
        kq = 1 + (f - 5) * (nq - 1)   ! k-index of the face in q-space

        if (interpolate) then
          associate(zf => z1(:,:,1))

            ! v⁻ at q-points ...................................................

            if (lobatto) then

              ! extraction from element q-points
              do c = 1, 3
                do j = 1, nq
                do i = 1, nq
                  vm_qf(i,j,c) = v_qe(i,j,kq,c)
                end do
                end do
              end do

            else

              ! interpolation from v-points to face q-points
              do c = 1, 3

                ! direction 1
                call LIBXSMM_DMMCall(xmm_set_f1, I_vq, v(:,:,kv,e,c), zf)

                ! direction 2
                call LIBXSMM_DMMCall(xmm_set_f2, zf, I_vq_t, vm_qf(:,:,c))

              end do

            end if

            ! v⁺ at q-points ...................................................

            ! interpolation from v-points to face q-points
            do c = 1, 3

              ! direction 1
              call LIBXSMM_DMMCall(xmm_set_f1, I_vq, vp(:,:,f,e,c), zf)

              ! direction 2
              call LIBXSMM_DMMCall(xmm_set_f2, zf, I_vq_t, vp_qf(:,:,c))

            end do

          end associate

        else ! no interpolation

          ! extract v⁻ and v⁺ ..................................................

          do c = 1, 3
            do j = 1, nv
            do i = 1, nv
              vm_qf(i,j,c) = v_qe(i,j,kv,c)
            end do
            end do
          end do

          do c = 1, 3
            do j = 1, nv
            do i = 1, nv
              vp_qf(i,j,c) = vp(i,j,f,e,c)
            end do
            end do
          end do

        end if

        ! integrand ............................................................

        select case(bc(f,e))

        case('D')

          ! boundary flux, assuming vb = {v} . . . . . . . . . . . . . . . . . .

          do j = 1, nq
          do i = 1, nq

            vb_1 = 0.5 * (vm_qf(i,j,1) + vp_qf(i,j,1))
            vb_2 = 0.5 * (vm_qf(i,j,2) + vp_qf(i,j,2))
            vb_3 = 0.5 * (vm_qf(i,j,3) + vp_qf(i,j,3))

            vb_n = n_q(i,j,f,e,1) * vb_1  &
                 + n_q(i,j,f,e,2) * vb_2  &
                 + n_q(i,j,f,e,3) * vb_3

            cf = -w_q(i) * w_q(j) * a_q(i,j,f,e) * vb_n

            ff(i,j,1) = cf * vb_1
            ff(i,j,2) = cf * vb_2
            ff(i,j,3) = cf * vb_3

          end do
          end do

        case default

          ! local Lax-Friedrichs flux  . . . . . . . . . . . . . . . . . . . . .

          do j = 1, nq
          do i = 1, nq

            vm_1 = vm_qf(i,j,1)
            vm_2 = vm_qf(i,j,2)
            vm_3 = vm_qf(i,j,3)
            vm_n = n_q(i,j,f,e,1) * vm_1  &
                 + n_q(i,j,f,e,2) * vm_2  &
                 + n_q(i,j,f,e,3) * vm_3

            vp_1 = vp_qf(i,j,1)
            vp_2 = vp_qf(i,j,2)
            vp_3 = vp_qf(i,j,3)
            vp_n = n_q(i,j,f,e,1) * vp_1  &
                 + n_q(i,j,f,e,2) * vp_2  &
                 + n_q(i,j,f,e,3) * vp_3

            cf = -w_q(i) * w_q(j) * a_q(i,j,f,e)
            vn = max(abs(vm_n), abs(vp_n))

            ff(i,j,1) = cf * (0.5 * (vm_n * vm_1 + vp_n * vp_1) + vn * (vm_1 - vp_1))
            ff(i,j,2) = cf * (0.5 * (vm_n * vm_2 + vp_n * vp_2) + vn * (vm_2 - vp_2))
            ff(i,j,3) = cf * (0.5 * (vm_n * vm_3 + vp_n * vp_3) + vn * (vm_3 - vp_3))

          end do
          end do

        end select

        ! face integral contributions ..........................................

        if (interpolate) then

          associate(zf => z1(:,:,1))
            do c = 1, 3

              ! integration over direction 2
              call LIBXSMM_DMMCall(xmm_set_f3, ff(:,:,c), I_vq, zf)

              ! integration over direction 1
              call LIBXSMM_DMMCall(xmm_add_f, I_vq_t, zf, F_c(:,:,kv,e,c))

            end do
          end associate

        else

          do c = 1, 3
            do j = 1, nv
            do i = 1, nv
              F_c(i,j,kv,e,c) = F_c(i,j,kv,e,c) + ff(i,j,c)
            end do
            end do
          end do

        end if

      end do Faces_5_and_6

    end do Elements


  end subroutine TPO_INS_Convection_D_XSMM_RDP

  !=============================================================================

#endif

end module TPO__INS_Convection__3D_D__XSMM_RDP
