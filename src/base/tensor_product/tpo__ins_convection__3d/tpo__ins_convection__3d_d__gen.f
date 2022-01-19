!> summary:  3d generic curvilinear element convection operator for INS (D)
!> author:   Jerome Michel, Jörg Stiller
!> date:     2021/12/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_INS_Convection_D_Gen &
               (nv, nq, ne, D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c)

  ! arguments ..................................................................

  integer, intent(in) :: nv !< number of velocity points per element direction
  integer, intent(in) :: nq !< number of quadrature points per direction
  integer, intent(in) :: ne !< number of elements

  real(RNP), intent(in) :: D_v (nv,nv) !< 1D diff operator on v-points
  real(RNP), intent(in) :: I_vq(nq,nv) !< 1D interpolation from to q-points
  real(RNP), intent(in) :: w_q (nq)    !< 1D quadrature weights

  real(RNP), intent(in)  :: Jd_q(nq,nq,nq,ne)     !< Jacobian determinant
  real(RNP), intent(in)  :: Ji_q(nq,nq,nq,ne,3,3) !< inverse Jacobian matrix
  real(RNP), intent(in)  :: a_q(nq,nq,6,ne)       !< face area coefficients
  real(RNP), intent(in)  :: n_q(nq,nq,6,ne)       !< face unit normal vectors
  real(RNP), intent(in)  :: v(nv,nv,nv,ne,3)      !< velocity
  real(RNP), intent(in)  :: vp(nv,nv,6,ne,3)      !< velocity v⁺ at faces
  real(RNP), intent(out) :: F_c(nv,nv,nv,ne,3)    !< convective flux integrals

  ! local variables ............................................................

  logical   :: interpolate      ! T/F if interpolation is/not required
  logical   :: lobatto          ! T/F if quadrature is/not of Lobatto type

  real(RNP) :: I_qv_t(nv,nq)    ! transposed interpolation operator
  real(RNP) :: D_vq(nq,nv)      ! interpolated differentiation operator
  real(RNP) :: M(nq,nq,nq)      ! mass matrix based on quadrature points

  real(RNP) :: v_qe(nq,nq,nq,3) ! velocity at element quadrature points
  real(RNP) :: vm_qf(nq,nq,3)   ! interior velocity v⁻ at face quadrature points
  real(RNP) :: vp_qf(nq,nq,3)   ! exterior velocity v⁺ at face quadrature points

  real(RNP) :: M_Ji_vv(nq,nq,nq,3,3) ! MJ⁻¹⋅vv at element quadrature points

  real(RNP) :: z1(nq,nv,nv)
  real(RNP) :: z2(nq,nq,nv)
  real(RNP) :: ff(nq,nq,3)

  real(RNP) :: v1v1, v1v2, v1v3, v2v2, v2v3, v3v3
  real(RNP) :: vm_1, vm_2, vm_3, vm_n
  real(RNP) :: vp_1, vp_2, vp_3, vp_n
  real(RNP) :: cf, vn, tmp

  integer   :: c, e, f, i, j, k, p
  integer   :: iv, iq, jv, jq, kv, kq

  !-----------------------------------------------------------------------------
  ! initialization

  lobatto = I_vq(1,1) == 1 .and. I_vq(nq,nv) == 1
  interpolate = nq /= nv  .or. .not. lobatto

  if (interpolate) then
    I_vq_t = transpose(I_vq)
    D_vq   = matmul(I_vq, D_v)
  else
    D_vq   = D_v
  end if

  !$omp do
  Elements: do e = 1, ne

    !---------------------------------------------------------------------------
    ! volume integral

    ! v at element quadrature points  ..........................................

    if (interpolate) then

      do c = 1, 3

        ! direction 1
        do k = 1, nv
        do j = 1, nv
        do i = 1, nq
          tmp = 0
          do p = 1, nv
            tmp = tmp + I_vq_t(p,i) * v(p,j,k,e,c)
          end do
          z1(i,j,k) = tmp
        end do
        end do
        end do

        ! direction 2
        do k = 1, nv
        do j = 1, nq
        do i = 1, nq
          tmp = 0
          do p = 1, nv
            tmp = tmp + I_vq_t(p,j) * z1(i,p,k)
          end do
          z2(i,j,k) = tmp
        end do
        end do
        end do

        ! direction 3
        do k = 1, nq
        do j = 1, nq
        do i = 1, nq
          tmp = 0
          do p = 1, nv
            tmp = tmp + I_vq_t(p,k) * z2(i,j,p)
          end do
          v_qe(i,j,k,c) = tmp
        end do
        end do
        end do

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

    ! mass matrix ..............................................................

    do k = 1, nq
    do j = 1, nq
    do i = 1, nq
      M(i,j,k) = w_q(k) * w_q(j) * w_q(i) * Jd_q(i,j,k,e)
    end do
    end do
    end do

    ! MJ⁻¹⋅vv ..................................................................

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

      ! note: reversed order to reuse workspace

      ! F_c = ∂₃ MJ⁻¹⋅vv (off-grid) ............................................

      do c = 1, 3

        ! projection in direction 3
        do k = 1, nv
        do j = 1, nq
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq(p,k) * M_Ji_vv(i,j,p,1,c)
          end do
          z2(i,j,k) = tmp
        end do
        end do
        end do

        ! projection in direction 2
        do k = 1, nv
        do j = 1, nv
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq(p,j) * z2(i,p,k)
          end do
          z1(i,j,k) = tmp
        end do
        end do
        end do

        ! differentiation in direction 1
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nq
            tmp = tmp + D_vq(p,i) * z1(p,j,k)
          end do
          F_c(i,j,k,e,c) = tmp
        end do
        end do
        end do

      end do

      ! F_c += ∂₂ MJ⁻¹⋅vv (off-grid) ...........................................

      do c = 1, 3

        ! projection in direction 3
        do k = 1, nv
        do j = 1, nq
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq_t(p,k) * M_Ji_vv(i,j,p,2,c)
          end do
          z2(i,j,k) = tmp
        end do
        end do
        end do

        ! differentiation in direction 2
        do k = 1, nv
        do j = 1, nv
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + D_vq(p,j) * z2(i,p,k)
          end do
          z1(i,j,k) = tmp
        end do
        end do
        end do

        ! projection in direction 1
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq(p,i) * z1(p,j,k)
          end do
          F_c(i,j,k,e,c) = F_c(i,j,k,e,c) + tmp
        end do
        end do
        end do

      end do

      ! F_c += ∂₁ MJ⁻¹⋅vv (off-grid) ...........................................

      do c = 1, 3

        ! differentiation in direction 3
        do k = 1, nv
        do j = 1, nq
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + D_vq(p,k) * M_Ji_vv(i,j,p,3,c)
          end do
          z2(i,j,k) = tmp
        end do
        end do
        end do

        ! projection in direction 2
        do k = 1, nv
        do j = 1, nv
        do i = 1, nq
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq(p,j) * z2(i,p,k)
          end do
          z1(i,j,k) = tmp
        end do
        end do
        end do

        ! projection in direction 1
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nq
            tmp = tmp + I_vq(p,i) * z1(p,j,k)
          end do
          F_c(i,j,k,e,c) = F_c(i,j,k,e,c) + tmp
        end do
        end do
        end do

      end do

    else

      ! standard divergence of MJ⁻¹⋅vv (on-grid) ...............................

      do  c = 1, 3

        ! direction 1
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nv
            tmp = tmp + D_vq(p,i) * M_Ji_vv(p,j,k,1,c)
          end do
          F_c(i,j,k,e,c) = tmp
        end do
        end do
        end do

        ! direction 2
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nv
            tmp = tmp + D_vq(p,j) * M_Ji_vv(i,p,k,2,c)
          end do
          F_c(i,j,k,e,c) = F_c(i,j,k,e,c) + tmp
        end do
        end do
        end do

        ! direction 3
        do k = 1, nv
        do j = 1, nv
        do i = 1, nv
          tmp = 0
          do p = 1, nv
            tmp = tmp + D_vq(p,k) * M_Ji_vv(i,j,p,3,c)
          end do
          F_c(i,j,k,e,c) = F_c(i,j,k,e,c) + tmp
        end do
        end do
        end do

      end do

    end if

    !---------------------------------------------------------------------------
    ! face integrals

    ! Γ₁ ∪ Γ₂ ..................................................................

    Faces_1_and_2: do f = 1, 2

      iv = 1 + (f - 1) * (nv - 1)   ! i-index of the face in v-space
      iq = 1 + (f - 1) * (nq - 1)   ! i-index of the face in q-space

      if (interpolate) then
        associate(zf => z1(:,:,1))

          ! v⁻ at q-points . . . . . . . . . . . . . . . . . . . . . . . . . . .

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
              do k = 1, nv
              do j = 1, nq
                tmp = 0
                do p = 1, nv
                  tmp = tmp + I_vq_t(p,j) * v(iv,p,k,e,c)
                end do
                zf(j,k) = tmp
              end do
              end do

              ! direction 3
              do k = 1, nq
              do j = 1, nq
                tmp = 0
                do p = 1, nv
                  tmp = tmp + I_vq_t(p,k) * zf(j,p)
                end do
                vm_qf(j,k,c) = tmp
              end do
              end do

            end do

          end if

          ! v⁺ at q-points . . . . . . . . . . . . . . . . . . . . . . . . . . .

          ! interpolation from v-points to face q-points
          do c = 1, 3

            ! direction 2
            do k = 1, nv
            do j = 1, nq
              tmp = 0
              do p = 1, nv
                tmp = tmp + I_vq_t(p,j) * vp(p,k,f,e,c)
              end do
              zf(j,k) = tmp
            end do
            end do

            ! direction 3
            do k = 1, nq
            do j = 1, nq
              tmp = 0
              do p = 1, nv
                tmp = tmp + I_vq_t(p,k) * zf(j,p)
              end do
              vp_qf(j,k,c) = tmp
            end do
            end do

          end do

        end associate

      else

        ! extract v⁻ and v⁺  . . . . . . . . . . . . . . . . . . . . . . . . . .

        do c = 1, 3
          do k = 1, nv
          do j = 1, nv
            vm_qf(j,k,c) = v_qe(iq,j,k,c)
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

      ! integrand  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      do k = 1, nq
      do j = 1, nq

        vm_1 = vm_qf(j,k,1)
        vm_2 = vm_qf(j,k,2)
        vm_3 = vm_qf(j,k,3)
        vm_n = n_q(j,k,f,e,1) * vm_1  &
             + n_q(j,k,f,e,2) * vm_2  &
             + n_q(j,k,f,e,2) * vm_3

        vp_1 = vp_qf(j,k,1)
        vp_2 = vp_qf(j,k,2)
        vp_3 = vp_qf(j,k,3)
        vp_n = n_q(j,k,f,e,1) * vp_1  &
             + n_q(j,k,f,e,2) * vp_2  &
             + n_q(j,k,f,e,2) * vp_3

        cf = -w_q(j) * w_q(k) * a_q(j,k,f,e)
        vn = max(abs(vm_n), abs(vp_n))

        ff(j,k,1) = cf * (0.5 * (vm_n * vm_1 + vp_n * vp_1) + vn * (vm_1 - vp_1))
        ff(j,k,2) = cf * (0.5 * (vm_n * vm_2 + vp_n * vp_2) + vn * (vm_2 - vp_2))
        ff(j,k,3) = cf * (0.5 * (vm_n * vm_3 + vp_n * vp_3) + vn * (vm_3 - vp_3))

      end do
      end do

      ! face integral contributions  . . . . . . . . . . . . . . . . . . . . . .

      if (interpolate) then

        associate(zf => z1(:,:,1))
          do c = 1, 3

            ! integration over direction 3
            do k = 1, nv
            do j = 1, nq
              tmp = 0
              do p = 1, nq
                tmp = tmp + I_vq(p,k) * ff(j,p,c)
              end do
              zf(j,k) = tmp
            end do
            end do

            ! integration over direction 2
            do k = 1, nv
            do j = 1, nv
              tmp = 0
              do p = 1, nq
                tmp = tmp + I_vq(p,k) * zf(p,k)
              end do
              F_c(iv,j,k,e,c) = F_c(iv,j,k,e,c) + tmp
            end do
            end do

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

    ! Γ₃ ∪ Γ₄ ..................................................................

    Faces_3_and_4: do f = 3, 4

      jv = 1 + (f - 3) * (nv - 1)   ! j-index of the face in v-space
      jq = 1 + (f - 3) * (nq - 1)   ! j-index of the face in q-space

      !### TBD ###

    end do Faces_3_and_4

    ! Γ₅ ∪ Γ₆ ..................................................................

    Faces_5_and_6: do f = 5, 6

      kv = 1 + (f - 5) * (nv - 1)   ! k-index of the face in v-space
      kq = 1 + (f - 5) * (nq - 1)   ! k-index of the face in q-space

      !### TBD ###

    end do Faces_5_and_6

  end do Elements

end subroutine TPO_INS_Convection_D_Gen
