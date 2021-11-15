!> summary:  3D DG diffusion operator: BC part of RHS with regular mesh and
!>           variable diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - An OpenMP race condition can appear if two or more faces of one element
!>     belong to the same boundary
!===============================================================================

submodule(DG__Diffusion_Operator__3D) MP_AddBC_RV
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Addition of BC to the RHS for regular mesh and variable ν
  !>
  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
  !> boundary values four each boundary. These values are applied to the
  !> right hand side `f` according boundary type specified in `this % bc`.

  module subroutine AddBC_RV(this, bv, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
    real(RNP), intent(inout) :: f(0:,0:,0:,:)

    ! local variables ..........................................................

    real(RNP), pointer, contiguous, save :: ub(:,:,:), qb(:,:,:)
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: Mfs
    real(RNP), dimension(0:this%eop%po, 6)             :: cd
    real(RNP), dimension(0:this%eop%po)                :: delta_0, delta_P

    real(RNP) :: ca(3), cx(3), mu(3)
    integer   :: b, e, i, j, k, l, m, n

    associate( mesh => this % sem % mesh &
             , nu   => this % nu_pv      &
             , Ms   => this % eop % w    &
             , Ds   => this % eop % D    &
             , P    => this % eop % po   )

      ! initialization .........................................................

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mfs(i,j) = Ms(i) * Ms(j)
      end do
      end do

      ! delta function
      delta_0 = [ ONE, (ZERO, i=1,P) ]
      delta_P = [ (ZERO, i=1,P), ONE ]

      ! penalties
      mu(1) = this % eop % PenaltyFactor(mesh % dx(1))
      mu(2) = this % eop % PenaltyFactor(mesh % dx(2))
      mu(3) = this % eop % PenaltyFactor(mesh % dx(3))

      ! metric coefficients
      ca(1) = mesh % dx(2) * mesh % dx(3) / 4
      ca(2) = mesh % dx(3) * mesh % dx(1) / 4
      ca(3) = mesh % dx(1) * mesh % dx(2) / 4
      cx(:) = 2 / mesh % dx(:)

      ! Dirichlet coefficients
      cd(:,1) = ca(1) * ( cx(1)*Ds(0,:) + 2*mu(1) * delta_0(:))
      cd(:,2) = ca(1) * (-cx(1)*Ds(P,:) + 2*mu(1) * delta_P(:))
      cd(:,3) = ca(2) * ( cx(2)*Ds(0,:) + 2*mu(2) * delta_0(:))
      cd(:,4) = ca(2) * (-cx(2)*Ds(P,:) + 2*mu(2) * delta_P(:))
      cd(:,5) = ca(3) * ( cx(3)*Ds(0,:) + 2*mu(3) * delta_0(:))
      cd(:,6) = ca(3) * (-cx(3)*Ds(P,:) + 2*mu(3) * delta_P(:))

      ! add boundary conditions ................................................

      Boundaries: do b = 1, mesh % n_bound

        Boundary_Faces: associate(face => mesh % boundary(b) % face)

          select case(this % bc(b))

          case('D')

            ! Dirichlet BC:  f = f - a/4 Mfs (n⋅ν⁻∇𝜑⁻ - 2μν 𝜑⁻) ub

            !$omp master
            ub(0:,0:,1:) => bv(b) % val(:,:,:,1)
            !$omp end master
            !$omp barrier

            !$omp do
            do l = 1, size(face)
              e = face(l) % mesh_element % id
              m = face(l) % mesh_element % face

              select case(m)

              case(1,2) ! direction 1
                n = (m-1) * P
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(i,m) * Mfs(j,k) * nu(n,j,k,e) * ub(j,k,l)
                end do
                end do
                end do

              case(3,4) ! direction 2
                n = (m-3) * P
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(j,m) * Mfs(i,k) * nu(i,n,k,e) * ub(i,k,l)
                end do
                end do
                end do

              case(5,6) ! direction 3
                n = (m-5) * P
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(k,m) * Mfs(i,j) * nu(i,j,n,e) * ub(i,j,l)
                end do
                end do
                end do

              end select
            end do

          case('N')

            ! Neumann BC:  f = f + a/4 Mfs 𝜑⁻ qb

            !$omp master
            qb(0:,0:,1:) => bv(b) % val(:,:,:,1)
            !$omp end master
            !$omp barrier

            !$omp do
            do l = 1, size(face)
              e = face(l) % mesh_element % id
              m = face(l) % mesh_element % face

              select case(m)

              case(1,2) ! direction 1
                i = (m-1) * P
                do k = 0, P
                do j = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + ca(1) * Mfs(j,k) * qb(j,k,l)
                end do
                end do

              case(3,4) ! direction 2
                j = (m-3) * P
                do k = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + ca(2) * Mfs(i,k) * qb(i,k,l)
                end do
                end do

              case(5,6) ! direction 3
                k = (m-5) * P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + ca(3) * Mfs(i,j) * qb(i,j,l)
                end do
                end do

              end select
            end do

          end select

        end associate Boundary_Faces
      end do Boundaries

    end associate

  end subroutine AddBC_RV

  !=============================================================================

end submodule MP_AddBC_RV
