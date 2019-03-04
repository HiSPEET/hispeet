!> summary:  ISP flow: computation of time derivative to given solution
!> author:   Joerg Stiller
!> date:     2018/05/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: computation of time derivative to given solution
!===============================================================================

module CART__ISP_Flow__Time_Derivative

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, TWO
  use Array_Assignments
  use TPO_sDDD
  use ISP_Flow_Problem
  use CART__TPO_Div
  use CART__TPO_Grad
  use CART__TPO_RotRot
  use CART__Mesh_Partition
  use CART__Weak_Gradient
  use CART__Weak_Divergence
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Convection

  implicit none
  private

  public :: TimeDerivative

contains

!-------------------------------------------------------------------------------
!> Computation of time derivative to given solution

subroutine TimeDerivative(problem, flow_op, t, u_c, u_d, p, nu, F, F_c, F_d, F_s)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(in)    :: flow_op !< flow operators
  real(RNP),            intent(in)    :: t       !< time
  real(RNP),  optional, intent(in)    :: u_c     !< variables for convection
  real(RNP),  optional, intent(in)    :: u_d     !< variables for diffusion
  real(RNP),  optional, intent(in)    :: p       !< pressure
  real(RNP),  optional, intent(in)    :: nu      !< variable diffusivity
  real(RNP),            intent(out)   :: F       !< ∂u/∂t
  real(RNP),  optional, intent(inout) :: F_c     !< convection part
  real(RNP),  optional, intent(inout) :: F_d     !< diffusion part
  real(RNP),  optional, intent(out)   :: F_s     !< source part

  dimension :: u_c  (:,:,:,:,:)
  dimension :: u_d  (:,:,:,:,:)
  dimension :: p    (:,:,:,:)
  dimension :: nu   (:,:,:,:,:)
  dimension :: F    (:,:,:,:,:)
  dimension :: F_c  (:,:,:,:,:)
  dimension :: F_d  (:,:,:,:,:)
  dimension :: F_s  (:,:,:,:,:)

  ! local variables ............................................................

  real(RNP), allocatable, save :: w(:,:,:,:,:)
  integer :: np, ne, nc

  associate( mesh => flow_op % mesh      &
           , Ms   => flow_op % eop_u % w &
           , Ds   => flow_op % eop_u % D &
           , Dd   => flow_op % eop_u % D )

    ! intialization ............................................................

    np = size(F,1)
    ne = size(F,4)
    nc = size(F,5)

    !$omp single
    allocate(w, mold=F)
    !$omp end single

    ! external sources .........................................................

    call problem % GetExternalSources(flow_op%x, t, F)
    if (present(F_s)) then
      call SetArray(F_s, F, multi=.true.)
    end if

    ! convection ...............................................................

    if (problem%stokes) then

      ! skip convection term
      if (present(F_c)) then
        call SetArray(F_c, ZERO, multi=.true.)
      end if

    else if (present(u_c)) then

      ! compute convection term
      call WeakConvectiveFlux(flow_op, u_c, div_F=w)
      call MergeArrays(ONE, F, -ONE, w, multi=.true.)
      if (present(F_c)) then
        call MergeArrays(ZERO, F_c, -ONE, w, multi=.true.)
      end if

    else if (present(F_c)) then

      ! use given convection term
      call MergeArrays(ONE, F, ONE, F_c, multi=.true.)

    end if

    ! diffusion ................................................................

    if (present(u_d)) then

      if (present(nu)) then
        call DiffTimeDeriv_Velocity_VI (mesh, Ms, Dd, nu, u_d, w, F, F_d)
        call DiffTimeDeriv_Scalars_VI  (mesh, Ms, Dd, nu, u_d, w, F, F_d)
      else
        associate(nu => problem % nu_ref)
          call DiffTimeDeriv_Velocity_CI (mesh, Ms, Dd, nu, u_d, w, F, F_d)
          call DiffTimeDeriv_Scalars_CI  (mesh, Ms, Dd, nu, u_d, w, F, F_d)
        end associate
      end if

      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,4), ZERO)
      end if

    else if (present(F_d)) then

      ! use given diffusion term
      call MergeArrays(ONE, F, ONE, F_d, multi=.true.)

    end if

    ! pressure .................................................................

    if (present(p)) then
      call WeakGradient(mesh, Ms, Ds, p, w)
      call MergeArrays(ONE, F, -ONE, w, multi=.true.)
    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(w)
    !$omp end master

  end associate

end subroutine TimeDerivative

!-------------------------------------------------------------------------------
!> Diffusive part of momentum equation with constant viscosity
!>
!> Computes
!>
!>     F(*,1:3) += ν∇·∇v - ν∇∇·v =F_d(*,1:3)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Velocity_CI(mesh, Ms, Dd, nu, v, w, F, F_d)
  class(MeshPartition), intent(in)    :: mesh           !< mesh partition
  real(RNP),            intent(in)    :: Ms (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu (:)         !< viscosity, ν = nu(1)
  real(RNP),            intent(in)    :: v  (:,:,:,:,:) !< velocity, v(*,1:3)
  real(RNP),            intent(inout) :: w  (:,:,:,:,:) !< workspace
  real(RNP),            intent(inout) :: F  (:,:,:,:,:) !< time derivative ∂v/∂t
  real(RNP),  optional, intent(out)   :: F_d(:,:,:,:,:) !< diff part of ∂v/∂t

  integer :: ne, np
  integer :: c

  np = size(Ms)
  ne = mesh%ne

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    ! div-grad part
    do c = 1, 3
      call TPO_Grad_Eval(np, ne, Dd, mesh%dx, v(:,:,:,:,c), g)
      call ScaleArray(g, nu(1), multi=.true.)
      call WeakDivergence(mesh, Ms, Dd, g, q)
      call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,c), q)
      end if
    end do

    ! divergence penalty
    call TPO_Div_Eval(np, ne, Dd, mesh%dx, v, q)
    call ScaleArray(q, -nu(1))
    call WeakGradient(mesh, Ms, Dd, q, g)
    call MergeArrays(ONE, F(:,:,:,:,1:3), ONE, g, multi=.true.)
    if (present(F_d)) then
      call MergeArrays(ONE, F_d(:,:,:,:,1:3), ONE, g, multi=.true.)
    end if

  end associate

end subroutine DiffTimeDeriv_Velocity_CI

!-------------------------------------------------------------------------------
!> Diffusive part of scalar transport with constant diffusivity
!>
!> Computes
!>
!>     F(*,5:nc) += ν∇·∇u = F_d(*,5:nc)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Scalars_CI(mesh, Ms, Dd, nu, u, w, F, F_d)
  class(MeshPartition), intent(in)    :: mesh           !< mesh partition
  real(RNP),            intent(in)    :: Ms (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu (:)         !< diffusivities
  real(RNP),            intent(in)    :: u  (:,:,:,:,:) !< scalars: u(*,5:nc)
  real(RNP),            intent(inout) :: w  (:,:,:,:,:) !< workspace
  real(RNP),            intent(inout) :: F  (:,:,:,:,:) !< time derivative ∂u/∂t
  real(RNP),  optional, intent(inout) :: F_d(:,:,:,:,:) !< diff part of ∂u/∂t

  procedure(TPO_Grad_Proc), pointer :: Gradient

  integer :: nc, ne, np
  integer :: c

  np = size(Ms)
  ne = size(u,4)
  nc = size(u,5)

  call TPO_Grad_Assign(np, Gradient)

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    do c = 5, nc
      call Gradient(np, ne, Dd, mesh%dx, u(:,:,:,:,c), w)
      call ScaleArray(g, nu(c), multi=.true.)
      call WeakDivergence(mesh, Ms, Dd, g, q)
      call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,c), q)
      end if
    end do

  end associate

end subroutine DiffTimeDeriv_Scalars_CI

!-------------------------------------------------------------------------------
!> Diffusive part of momentum equation with variable viscosity
!>
!> Computes
!>
!>     F(*,1:3) += ∇·[ν∇v + ∇·ν(∇v)ᵀ] - 2ν∇∇·v = F_d(*,1:3)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Velocity_VI(mesh, Ms, Dd, nu, v, w, F, F_d)
  class(MeshPartition), intent(in)    :: mesh           !< mesh partition
  real(RNP),            intent(in)    :: Ms (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu (:,:,:,:,:) !< viscosity, ν = nu(*,1)
  real(RNP),            intent(in)    :: v  (:,:,:,:,:) !< velocity, v(*,1:3)
  real(RNP),            intent(inout) :: w  (:,:,:,:,:) !< workspace
  real(RNP),            intent(inout) :: F  (:,:,:,:,:) !< time derivative ∂v/∂t
  real(RNP),  optional, intent(out)   :: F_d(:,:,:,:,:) !< diff part of ∂v/∂t

  real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:)
  integer :: ne, np
  integer :: c, d, e, i, j, k

  np = size(Ms)
  ne = mesh%ne

  !$omp single
  allocate(grad_v(np,np,np,ne,3,3))
  !$omp end single

  call TPO_Grad_Eval(np, 3*ne, Dd, mesh%dx, v(:,:,:,:,1:3), grad_v)

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    ! diffusion part ...........................................................

    do c = 1, 3

      do d = 1, 3
        !$omp do
        do e = 1, ne
          do k = 1, np
          do j = 1, np
          do i = 1, np
            g(i,j,k,e,d) = nu(i,j,k,e,1) * ( grad_v(i,j,k,e,d,c) &
                                           + grad_v(i,j,k,e,c,d) )
          end do
          end do
          end do
        end do
      end do

      call WeakDivergence(mesh, Ms, Dd, g, q)

      call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,c), q)
      end if

    end do

    ! div penalty ..............................................................

    !$omp do
    do e = 1, ne
      do k = 1, np
      do j = 1, np
      do i = 1, np
        q(i,j,k,e) = grad_v(i,j,k,e,1,1) &
                   + grad_v(i,j,k,e,2,2) &
                   + grad_v(i,j,k,e,3,3)
      end do
      end do
      end do
    end do

    call WeakGradient(mesh, Ms, Dd, q, g)

    do d = 1, 3
      !$omp do
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np
          g(i,j,k,e,d) = nu(i,j,k,e,1) * g(i,j,k,e,d)
        end do
        end do
        end do
      end do
    end do

    call MergeArrays(ONE, F, -TWO, g, multi = .true.)
    if (present(F_d)) then
      call SetArray(F_d, g, multi = .true.)
    end if

  end associate

  !$omp barrier
  !$omp master
  deallocate(grad_v)
  !$omp end master

end subroutine DiffTimeDeriv_Velocity_VI

!-------------------------------------------------------------------------------
!> Diffusive part of scalar transport with variable diffusivity
!>
!> Computes
!>
!>     F(*,5:nc) += ∇·ν∇u = F_d(*,5:nc)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Scalars_VI(mesh, Ms, Dd, nu, u, w, F, F_d)
  class(MeshPartition), intent(in)    :: mesh           !< mesh partition
  real(RNP),            intent(in)    :: Ms (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu (:,:,:,:,:) !< diffusivities
  real(RNP),            intent(in)    :: u  (:,:,:,:,:) !< scalars: u(*,5:nc)
  real(RNP),            intent(inout) :: w  (:,:,:,:,:) !< workspace
  real(RNP),            intent(inout) :: F  (:,:,:,:,:) !< time derivative ∂u/∂t
  real(RNP),  optional, intent(inout) :: F_d(:,:,:,:,:) !< diff part of ∂u/∂t

  procedure(TPO_Grad_Proc), pointer :: Gradient

  integer :: nc, ne, np
  integer :: c, d, e, i, j, k

  np = size(Ms)
  ne = size(u,4)
  nc = size(u,5)

  call TPO_Grad_Assign(np, Gradient)

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    do c = 5, nc

      call Gradient(np, ne, Dd, mesh%dx, u(:,:,:,:,c), g)

      do d = 1, 3
        !$omp do
        do e = 1, ne
          do k = 1, np
          do j = 1, np
          do i = 1, np
            g(i,j,k,e,d) = nu(i,j,k,e,1) * g(i,j,k,e,d)
          end do
          end do
          end do
        end do
      end do

      call WeakDivergence(mesh, Ms, Dd, g, q)

      call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,c), q)
      end if

    end do

  end associate

end subroutine DiffTimeDeriv_Scalars_VI

!===============================================================================

end module CART__ISP_Flow__Time_Derivative