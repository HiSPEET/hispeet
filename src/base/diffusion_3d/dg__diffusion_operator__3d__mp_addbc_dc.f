!> summary:  3D DG diffusion operator: BC part of RHS with deformed mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/10/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - An OpenMP race condition can appear if two or more faces of one element
!>     belong to the same boundary
!===============================================================================

submodule(DG__Diffusion_Operator__3D) MP_AddBC_DC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Addition of BC to the RHS for deformed mesh and constant ν
  !>
  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
  !> boundary values four each boundary. These values are applied to the
  !> right hand side `f` according boundary type specified in `this % bc`.

  module subroutine AddBC_DC(this, bv, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
    real(RNP), intent(inout) :: f(0:,0:,0:,:)

    ! local variables ..........................................................

    real(RNP), pointer, contiguous, save :: ub(:,:,:), qb(:,:,:)
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po, 3) :: Cu
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po)    :: Mfs
    real(RNP) :: mu_nu, tmp
    integer   :: b, e, i, j, k, l, m, n

    associate( mesh => this % sem % mesh           &
             , a    => this % sem % metrics % a    &
             , Ji_n => this % sem % metrics % Ji_n &
             , nu   => this % nu_pc                &
             , Ms   => this % eop % w              &
             , Ds   => this % eop % D              &
             , P    => this % eop % po             )

      ! auxiliaries ............................................................

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mfs(i,j) = Ms(i) * Ms(j)
      end do
      end do

      ! add boundary conditions ................................................

      Boundaries: do b = 1, mesh % n_bound

        Boundary_Faces: associate(face => mesh % boundary(b) % face)

          select case(this % bc(b))

          case('D')

            ! Dirichlet BC . . . . . . . . . . . . . . . . . . . . . . . . . . .

            !$omp master
            ub(0:,0:,1:) => bv(b) % val(:,:,:,1)
            !$omp end master
            !$omp barrier

            !$omp do
            do l = 1, size(face)
              e = face(l) % mesh_element % id
              m = face(l) % mesh_element % face

              mu_nu = this % eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu

              select case(m)

              case(1,2)

                i = (m-1) * P

                do k = 0, P
                do j = 0, P
                  tmp = -Mfs(j,k) * nu * a(j,k,m,e) * ub(j,k,l)
                  Cu(j,k,1) = tmp * Ji_n(j,k,m,e,1)
                  Cu(j,k,2) = tmp * Ji_n(j,k,m,e,2)
                  Cu(j,k,3) = tmp * Ji_n(j,k,m,e,3)
                end do
                end do

                do k = 0, P
                do j = 0, P

                  do n = 0, P
                    f(n,j,k,e) = f(n,j,k,e) + Ds(i,n) * Cu(j,k,1)
                  end do

                  tmp = Mfs(j,k) * a(j,k,m,e) * mu_nu * 2 * ub(j,k,l)
                  do n = 0, P
                    tmp = tmp + Ds(n,j) * Cu(n,k,2) + Ds(n,k) * Cu(j,n,3)
                  end do
                  f(i,j,k,e) = f(i,j,k,e) + tmp

                end do
                end do

              case(3,4)

                j = (m-3) * P

                do k = 0, P
                do i = 0, P
                  tmp = -Mfs(i,k) * nu * a(i,k,m,e) * ub(i,k,l)
                  Cu(i,k,1) = tmp * Ji_n(i,k,m,e,1)
                  Cu(i,k,2) = tmp * Ji_n(i,k,m,e,2)
                  Cu(i,k,3) = tmp * Ji_n(i,k,m,e,3)
                end do
                end do

                do k = 0, P
                do i = 0, P

                  do n = 0, P
                    f(i,n,k,e) = f(i,n,k,e) + Ds(j,n) * Cu(i,k,2)
                  end do

                  tmp = Mfs(i,k) * a(i,k,m,e) * mu_nu * 2 * ub(i,k,l)
                  do n = 0, P
                    tmp = tmp + Ds(n,i) * Cu(n,k,1) + Ds(n,k) * Cu(i,n,3)
                  end do
                  f(i,j,k,e) = f(i,j,k,e) + tmp

                end do
                end do

              case(5,6)

                k = (m-5) * P

                do j = 0, P
                do i = 0, P
                  tmp = -Mfs(i,j) * nu * a(i,j,m,e) * ub(i,j,l)
                  Cu(i,j,1) = tmp * Ji_n(i,j,m,e,1)
                  Cu(i,j,2) = tmp * Ji_n(i,j,m,e,2)
                  Cu(i,j,3) = tmp * Ji_n(i,j,m,e,3)
                end do
                end do

                do j = 0, P
                do i = 0, P

                  do n = 0, P
                    f(i,j,n,e) = f(i,j,n,e) + Ds(k,n) * Cu(i,j,3)
                  end do

                  tmp = Mfs(i,j) * a(i,j,m,e) * mu_nu * 2 * ub(i,j,l)
                  do n = 0, P
                    tmp = tmp + Ds(n,i) * Cu(n,j,1) + Ds(n,j) * Cu(i,n,2)
                  end do
                  f(i,j,k,e) = f(i,j,k,e) + tmp

                end do
                end do

              end select
            end do
            !$omp end do nowait

          case('N')

            ! Neumann BC . . . . . . . . . . . . . . . . . . . . . . . . . . . .

            !$omp master
            qb(0:,0:,1:) => bv(b) % val(:,:,:,1)
            !$omp end master
            !$omp barrier

            !$omp do
            do l = 1, size(face)
              e = face(l) % mesh_element % id
              m = face(l) % mesh_element % face

              select case(m)

              case(1,2)
                i = (m-1) * P
                do k = 0, P
                do j = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + Mfs(j,k) * a(j,k,m,e) * qb(j,k,l)
                end do
                end do

              case(3,4)
                j = (m-3) * P
                do k = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + Mfs(i,k) * a(i,k,m,e) * qb(i,k,l)
                end do
                end do

              case(5,6)
                k = (m-5) * P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + Mfs(i,j) * a(i,j,m,e) * qb(i,j,l)
                end do
                end do

              end select
            end do
            !$omp end do nowait

          end select

        end associate Boundary_Faces
      end do Boundaries
      !$omp barrier

    end associate

  end subroutine AddBC_DC

  !=============================================================================

end submodule MP_AddBC_DC
