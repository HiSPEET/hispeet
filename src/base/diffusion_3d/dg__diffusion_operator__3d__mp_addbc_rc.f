!> summary:  3D DG diffusion operator: BC part of RHS with regular mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - An OpenMP race condition can appear if two or more faces of one element
!>     belong to the same boundary
!===============================================================================

submodule(DG__Diffusion_Operator__3D) MP_AddBC_RC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Addition of BC to the RHS for regular mesh and constant ν
  !>
  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
  !> boundary values four each boundary. These values are applied to the
  !> right hand side `f` according boundary type specified in `this % bc`.

  module subroutine AddBC_RC(this, bv, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
    real(RNP), intent(inout) :: f(:,:,:,:)

    ! local variables ..........................................................

    real(RNP), pointer, contiguous, save :: ub(:,:,:), qb(:,:,:)
    real(RNP), dimension(0:this%eop%po, 0:this%eop%po) :: As, Bs, Mfs
    real(RNP), dimension(0:this%eop%po, 6)             :: cd
    real(RNP), dimension(0:this%eop%po)                :: delta_0, delta_P

    real(RNP) :: ca(3), cx(3), mu(3)
    integer   :: b, e, i, j, k, l, m, P

    associate( mesh => this % sem % mesh &
             , nu_p => this % nu_pc      &
             , nu_s => this % nu_sc      &
             , eop  => this % eop        )

      ! initialization .........................................................

      P = eop % po

      ! computation of 1D standard diffusion and standard flux operator
      if (nu_s > 0) then
        call eop % Get_SVV_StandardStiffnessMatrix(As)
        call eop % Get_SVV_StandardDiffMatrix(Bs)
        As = nu_s * As
        Bs = nu_s * Bs
      else
        As = ZERO
        Bs = ZERO
      end if

      As = As + nu_p * eop%L
      Bs = Bs + nu_p * eop%D

      ! face quadrature weights
      do j = 0, P
      do i = 0, P
        Mfs(i,j) = eop % w(i) * eop % w(j)
      end do
      end do

      ! delta function
      delta_0 = [ ONE, (ZERO, i=1,P) ]
      delta_P = [ (ZERO, i=1,P), ONE ]

      ! penalties
      mu(1) = eop % PenaltyFactor(mesh % dx(1))
      mu(2) = eop % PenaltyFactor(mesh % dx(2))
      mu(3) = eop % PenaltyFactor(mesh % dx(3))

      ! metric coefficients
      ca(1) = mesh % dx(2) * mesh % dx(3) / 4
      ca(2) = mesh % dx(3) * mesh % dx(1) / 4
      ca(3) = mesh % dx(1) * mesh % dx(2) / 4
      cx(:) = 2 / mesh % dx(:)

      ! Dirichlet coefficients
      cd(:,1) = ca(1) * ( cx(1)*Bs(0,:) + 2*mu(1) * (nu_p + nu_s) * delta_0(:))
      cd(:,2) = ca(1) * (-cx(1)*Bs(P,:) + 2*mu(1) * (nu_p + nu_s) * delta_P(:))
      cd(:,3) = ca(2) * ( cx(2)*Bs(0,:) + 2*mu(2) * (nu_p + nu_s) * delta_0(:))
      cd(:,4) = ca(2) * (-cx(2)*Bs(P,:) + 2*mu(2) * (nu_p + nu_s) * delta_P(:))
      cd(:,5) = ca(3) * ( cx(3)*Bs(0,:) + 2*mu(3) * (nu_p + nu_s) * delta_0(:))
      cd(:,6) = ca(3) * (-cx(3)*Bs(P,:) + 2*mu(3) * (nu_p + nu_s) * delta_P(:))

      ! add boundary conditions ................................................

      Boundaries: do b = 1, mesh % n_bound

        Boundary_Faces: associate(face => mesh % boundary(b) % face)

          select case(this % bc(b))

          case('D')

            ! Dirichlet BC:  f = f - a/4 Mfs (n⋅(ν+νˢQ)∇𝜑⁻ - 2μ(ν+νˢ)𝜑⁻) ub

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
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(i,m) * Mfs(j,k) * ub(j,k,l)
                end do
                end do
                end do

              case(3,4) ! direction 2
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(j,m) * Mfs(i,k) * ub(i,k,l)
                end do
                end do
                end do

              case(5,6) ! direction 3
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(k,m) * Mfs(i,j) * ub(i,j,l)
                end do
                end do
                end do

              end select
            end do
            !$omp end do nowait

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
            !$omp end do nowait

          end select

        end associate Boundary_Faces
      end do Boundaries
      !$omp barrier

    end associate

  end subroutine AddBC_RC

  !=============================================================================

end submodule MP_AddBC_RC
