!> summary:  Incompressible Navier-Stokes viscous stress vector on boundary (DV)
!> author:   Joerg Stiller
!> date:     2024/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetViscousBoundaryStress_V

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Viscous stress vector on a boundary (V)

  module subroutine GetViscousBoundaryStress_V &
      (this, b, mu, nu, v, sb, bv, xout, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    integer, intent(in) :: b
    !< boundary ID

    real(RNP), contiguous, intent(in) :: mu(0:,0:,0:,:)
    !< kinematic bulk viscosity μ

    real(RNP), contiguous, intent(in) :: nu(0:,0:,0:,:)
    !< kinematic shear viscosity ν

    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)
    !< velocity

    real(RNP), contiguous, intent(out) :: sb(0:,0:,:,:)
    !< viscous stress vector on boundary `b` (0:po,0:po,1:n_face,1:3)

    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values at final time t
    !!   - Γᴰ :  vᵇ    in components 1:3, currently unused
    !!   - Γᴼ :  τ_nn  in component    4, including penalty term
    !!
    !! ignored at extrapolated boundaries

    logical, optional, intent(in) :: xout
    !< ignore boundary values and extrapolate `∇v` [F]

    integer, optional, intent(in) :: form
    !< switch form of the stress tensor `τ`
    !!   - `0` :  `ν[∇v + (∇v)ᵀ] + (μ-2/3ν)(∇⋅v)I`   (default)
    !!   - `1` :  `ν[∇v + (∇v)ᵀ] -       ν (∇⋅v)I`   (diffusion)
    !!   - `2` :  `ν[∇v + (∇v)ᵀ] -      2ν (∇⋅v)I`   (rotation)
    !!   - `3` :  `                      ν (∇⋅v)I`   (divergence penalty)

    ! local variables ..........................................................

    real(RNP), allocatable :: grad_v(:,:,:,:)
    real(RNP), allocatable :: chi_f(:,:), mu_f(:,:), nu_f(:,:)
    real(RNP) :: c_mu, c_nu, d_nu
    logical   :: outflow_bc, outflow_bv
    integer   :: e, f, m

    ! initialization ...........................................................

    if (present(xout)) then
      outflow_bc = this % problem % bc_v(b) == 'O' .and. .not. xout
    else
      outflow_bc = this % problem % bc_v(b) == 'O'
    end if
    outflow_bv = outflow_bc .and. present(bv)

    if (present(form)) then
      select case(form)
      case(1)
        c_mu =  0
        c_nu = -1
        d_nu =  1
      case(2)
        c_mu =  0
        c_nu = -2
        d_nu =  1
      case(3)
        c_mu =  0
        c_nu =  1
        d_nu =  0
      case default
        c_mu =  1
        c_nu = -2 * THIRD
        d_nu =  1
      end select
    else
      c_mu =  1
      c_nu = -2 * THIRD
      d_nu =  1
    end if

    associate( boundary => this % sem_u % mesh % boundary(b) &
             , Ji       => this % sem_u % metrics % Ji       &
             , n        => this % sem_u % metrics % n        &
             , po       => this % eop_u % po                 &
             , Ds       => this % eop_u % D                  )

      ! workspace ..............................................................

      allocate(grad_v(0:po,0:po,3,3), chi_f(0:po,0:po))
      allocate(mu_f, nu_f, mold = chi_f)

      !$omp do
      do f = 1, boundary % n_face
        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

        call GetVelocityGradient(e, m, po, Ds, Ji, v, grad_v)

        if (outflow_bc) then
          call ApplyOutflowBC(e, m, po, n, grad_v)
        end if

        select case(m)
        case(1)
          mu_f = mu( 0,:,:,e)
          nu_f = nu( 0,:,:,e)
        case(2)
          mu_f = mu(po,:,:,e)
          nu_f = nu(po,:,:,e)
        case(3)
          mu_f = mu(:, 0,:,e)
          nu_f = nu(:, 0,:,e)
        case(4)
          mu_f = mu(:,po,:,e)
          nu_f = nu(:,po,:,e)
        case(5)
          mu_f = mu(:,:, 0,e)
          nu_f = nu(:,:, 0,e)
        case(6)
          mu_f = mu(:,:,po,e)
          nu_f = nu(:,:,po,e)
        end select
        chi_f = c_mu * nu_f + c_nu * nu_f
        nu_f  = d_nu * nu_f

        call GetViscousStressVector(f, e, m, po, n, chi_f, nu_f, grad_v, sb)

        if (outflow_bv) then
          associate(tau_nn => bv(b) % val(:,:,f,4))
            sb(:,:,f,1) = sb(:,:,f,1) + n(:,:,m,e,1) * tau_nn
            sb(:,:,f,2) = sb(:,:,f,2) + n(:,:,m,e,2) * tau_nn
            sb(:,:,f,3) = sb(:,:,f,3) + n(:,:,m,e,3) * tau_nn
          end associate
        end if

      end do

    end associate

  end subroutine GetViscousBoundaryStress_V

  !-----------------------------------------------------------------------------
  !> Computation of the velocity gradient on a given element face

  pure subroutine GetVelocityGradient(e, m, po, Ds, Ji, v, grad_v)

    integer, intent(in) :: e
    !< element ID
    integer, intent(in) :: m
    !< element face
    integer, intent(in) :: po
    !< polynomial order
    real(RNP), contiguous, intent(in) :: Ds(0:,0:)
    !< standard diff matrix
    real(RNP), contiguous, intent(in) :: Ji(0:,0:,0:,:,:,:)
    !< inverse element Jacobian matrix
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)
    !< velocity
    real(RNP), contiguous, intent(out) :: grad_v(0:,0:,:,:)
    !< velocity gradient

    real(RNP) :: w(0:po,0:po,3)
    integer   :: c, i, j, k, p

    select case(m)

    case(1:2)

      i = (m-1) * po

      do c = 1, 3

        ! standard derivatives of component c
        w = 0
        do k = 0, po
        do j = 0, po
          do p = 0, po
            w(j,k,1) = w(j,k,1) + Ds(i,p) * v(p,j,k,e,c)
            w(j,k,2) = w(j,k,2) + Ds(j,p) * v(i,p,k,e,c)
            w(j,k,3) = w(j,k,3) + Ds(k,p) * v(i,j,p,e,c)
          end do
        end do
        end do

        ! apply inverse Jacobian matrix
        do k = 0, po
        do j = 0, po

          grad_v(j,k,1,c) = w(j,k,1) * Ji(i,j,k,e,1,1) &
                          + w(j,k,2) * Ji(i,j,k,e,2,1) &
                          + w(j,k,3) * Ji(i,j,k,e,3,1)

          grad_v(j,k,2,c) = w(j,k,1) * Ji(i,j,k,e,1,2) &
                          + w(j,k,2) * Ji(i,j,k,e,2,2) &
                          + w(j,k,3) * Ji(i,j,k,e,3,2)

          grad_v(j,k,3,c) = w(j,k,1) * Ji(i,j,k,e,1,3) &
                          + w(j,k,2) * Ji(i,j,k,e,2,3) &
                          + w(j,k,3) * Ji(i,j,k,e,3,3)

        end do
        end do

      end do

    case(3:4)

      j = (m-3) * po

      do c = 1, 3

        ! standard derivatives of component c
        w = 0
        do k = 0, po
        do i = 0, po
          do p = 0, po
            w(i,k,1) = w(i,k,1) + Ds(i,p) * v(p,j,k,e,c)
            w(i,k,2) = w(i,k,2) + Ds(j,p) * v(i,p,k,e,c)
            w(j,k,3) = w(j,k,3) + Ds(k,p) * v(i,j,p,e,c)
          end do
        end do
        end do

        ! apply inverse Jacobian matrix
        do k = 0, po
        do i = 0, po

          grad_v(i,k,1,c) = w(i,k,1) * Ji(i,j,k,e,1,1) &
                          + w(i,k,2) * Ji(i,j,k,e,2,1) &
                          + w(i,k,3) * Ji(i,j,k,e,3,1)

          grad_v(i,k,2,c) = w(i,k,1) * Ji(i,j,k,e,1,2) &
                          + w(i,k,2) * Ji(i,j,k,e,2,2) &
                          + w(i,k,3) * Ji(i,j,k,e,3,2)

          grad_v(i,k,3,c) = w(i,k,1) * Ji(i,j,k,e,1,3) &
                          + w(i,k,2) * Ji(i,j,k,e,2,3) &
                          + w(i,k,3) * Ji(i,j,k,e,3,3)

        end do
        end do

      end do

    case(5:6)

      k = (m-5) * po

      do c = 1, 3

        ! standard derivatives of component c
        w = 0
        do j = 0, po
        do i = 0, po
          do p = 0, po
            w(i,j,1) = w(i,j,1) + Ds(i,p) * v(p,j,k,e,c)
            w(j,k,2) = w(j,k,2) + Ds(j,p) * v(i,p,k,e,c)
            w(j,k,3) = w(j,k,3) + Ds(k,p) * v(i,j,p,e,c)
          end do
        end do
        end do

        ! apply inverse Jacobian matrix
        do j = 0, po
        do i = 0, po

          grad_v(i,j,1,c) = w(i,j,1) * Ji(i,j,k,e,1,1) &
                          + w(i,j,2) * Ji(i,j,k,e,2,1) &
                          + w(i,j,3) * Ji(i,j,k,e,3,1)

          grad_v(i,j,2,c) = w(i,j,1) * Ji(i,j,k,e,1,2) &
                          + w(i,j,2) * Ji(i,j,k,e,2,2) &
                          + w(i,j,3) * Ji(i,j,k,e,3,2)

          grad_v(i,j,3,c) = w(i,j,1) * Ji(i,j,k,e,1,3) &
                          + w(i,j,2) * Ji(i,j,k,e,2,3) &
                          + w(i,j,3) * Ji(i,j,k,e,3,3)

        end do
        end do

      end do

    end select

  end subroutine GetVelocityGradient

  !-----------------------------------------------------------------------------
  !> Computation of the velocity gradient on a given element face

  pure subroutine ApplyOutflowBC(e, m, po, n, grad_v)

    integer, intent(in) :: e
    !< element ID
    integer, intent(in) :: m
    !< element face
    integer, intent(in) :: po
    !< polynomial order
    real(RNP), contiguous, intent(in) :: n(0:,0:,:,:,:)
    !< element face normal vector
    real(RNP), contiguous, intent(inout) :: grad_v(0:,0:,:,:)
    !< velocity gradient, ∇v ← ∇v - nn⋅∇v

    real(RNP) :: dn_v(0:po,0:po)
    integer   :: i

    ! remove normal derivative: ∇v = ∇v - n(n⋅∇v)
    do i = 1, 3

      ! ∂vᵢ/∂n = n⋅∇vᵢ
      dn_v = n(:,:,m,e,1) * grad_v(:,:,1,i) &
           + n(:,:,m,e,2) * grad_v(:,:,2,i) &
           + n(:,:,m,e,3) * grad_v(:,:,3,i)

      ! ∇vᵢ = ∇vᵢ - n ∂vᵢ/∂n
      grad_v(:,:,1,i) = grad_v(:,:,1,i) - n(:,:,m,e,1) * dn_v
      grad_v(:,:,2,i) = grad_v(:,:,2,i) - n(:,:,m,e,2) * dn_v
      grad_v(:,:,3,i) = grad_v(:,:,3,i) - n(:,:,m,e,3) * dn_v

    end do

  end subroutine ApplyOutflowBC

  !-----------------------------------------------------------------------------
  !> Computation of the stress vector on a given element face

  pure subroutine GetViscousStressVector(f, e, m, po, n, chi, nu, grad_v, sb)

    integer, intent(in) :: f
    !< boundary face ID
    integer, intent(in) :: e
    !< element ID
    integer, intent(in) :: m
    !< element face
    integer, intent(in) :: po
    !< polynomial order
    real(RNP), contiguous, intent(in) :: n(0:,0:,:,:,:)
    !< element face normal vector
    real(RNP), contiguous, intent(in) :: chi(0:,0:)
    !< bulk viscosity, χ = (ζ - 2/3 η) / ρ
    real(RNP), contiguous, intent(in) :: nu(0:,0:)
    !< shear viscosity, ν = η/ρ
    real(RNP), contiguous, intent(in) :: grad_v(0:,0:,:,:)
    !< velocity gradient
    real(RNP), contiguous, intent(out) :: sb(0:,0:,:,:)
    !< boundary face stress vector

    real(RNP) :: tau(3,3), div_v
    integer   :: c, i, j

    do j = 0, po
    do i = 0, po

      ! divergence
      div_v = grad_v(i,j,1,1) + grad_v(i,j,2,2) + grad_v(i,j,3,3)

      ! stress tensor
      tau(1,1) = nu(i,j) * 2 * grad_v(i,j,1,1) + chi(i,j) * div_v
      tau(2,2) = nu(i,j) * 2 * grad_v(i,j,2,2) + chi(i,j) * div_v
      tau(3,3) = nu(i,j) * 2 * grad_v(i,j,3,3) + chi(i,j) * div_v
      tau(1,2) = nu(i,j) * (grad_v(i,j,1,2) + grad_v(i,j,2,1))
      tau(1,3) = nu(i,j) * (grad_v(i,j,1,3) + grad_v(i,j,3,1))
      tau(2,3) = nu(i,j) * (grad_v(i,j,2,3) + grad_v(i,j,3,2))

      ! stress vector
      do c = 1, 3
        sb(i,j,f,c) = n(i,j,m,e,1) * tau(1,c) &
                    + n(i,j,m,e,2) * tau(2,c) &
                    + n(i,j,m,e,3) * tau(3,c)
      end do

    end do
    end do

  end subroutine GetViscousStressVector

  !=============================================================================

end submodule MP_GetViscousBoundaryStress_V
