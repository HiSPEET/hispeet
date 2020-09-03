!> summary:  ISP flow: computation of time derivative to given solution
!> author:   Joerg Stiller
!> date:     2018/05/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Computes the time derivative to a given solution and/or contributions to it.
!>
!>       F = ∂u/∂t = F_c + F_d + F_p + F_s
!>
!> where
!>
!>      F_c = -∇·(v u_c)
!>
!> with `v = u_c(*,1:3)` if `u_c` is given and `F_c = 0` else,
!>
!>      F_d = F_d1 + F_d2 + F_d3
!>
!>      F_d1(*, : ) = ∇·(ν (∇v) )
!>      F_d2(*,1:3) = ∇·(ν (∇v)ᵀ)
!>      F_d3(*,1:3) = χ∇(ν (∇·v))
!>
!> with `v = u_d(*,1:3)` if `u_d` is given and `F_d = F_d1 = Fd2 = F_d3 = 0`
!> else,
!>
!>      F_p = -∇p
!>
!> if `p` is present and `F_p = 0` else, and
!>
!>      F_s = f(x,t)
!>
!> are the solution-independent sources.
!>
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
  use CART__Boundary_Variable
  use CART__DG_Weak_Gradient
  use CART__DG_Weak_Divergence
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Convection

  implicit none
  private

  public :: TimeDerivative

contains

!-------------------------------------------------------------------------------
!> Computation of time derivative to given solution
!>
!> @note
!>   * if absent, the bulk diffusion parameter χ ist taken from `flow_op`
!> @endnote

subroutine TimeDerivative( problem, flow_op, t, u_c, u_d, p, nu, chi,      &
                           F, F_c, F_d, F_d1, F_d2, F_d3, F_p, F_s, source )

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(in)    :: flow_op !< flow operators
  real(RNP),            intent(in)    :: t       !< time
  real(RNP),  optional, intent(in)    :: u_c     !< variables for convection
  real(RNP),  optional, intent(in)    :: u_d     !< variables for diffusion
  real(RNP),  optional, intent(in)    :: p       !< pressure
  real(RNP),  optional, intent(in)    :: nu      !< variable diffusivity
  real(RNP),  optional, intent(in)    :: chi     !< bulk diffusion parameter χ
  real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t
  real(RNP),  optional, intent(inout) :: F_c     !< convection part
  real(RNP),  optional, intent(inout) :: F_d     !< diffusion part, complete
  real(RNP),  optional, intent(out)   :: F_d1    !< diffusion, only ∇·ν∇u
  real(RNP),  optional, intent(out)   :: F_d2    !< diffusion, only ∇·ν(∇v)ᵀ
  real(RNP),  optional, intent(out)   :: F_d3    !< diffusion, only χ∇ν(∇·v)
  real(RNP),  optional, intent(out)   :: F_p     !< pressure part
  real(RNP),  optional, intent(out)   :: F_s     !< source part
  logical,    optional, intent(in)    :: source  !< include sources in F [T]

  dimension :: u_c  (:,:,:,:,:)
  dimension :: u_d  (:,:,:,:,:)
  dimension :: p    (:,:,:,:)
  dimension :: nu   (:,:,:,:,:)
  dimension :: F    (:,:,:,:,:)
  dimension :: F_c  (:,:,:,:,:)
  dimension :: F_d  (:,:,:,:,:)
  dimension :: F_d1 (:,:,:,:,:)
  dimension :: F_d2 (:,:,:,:,:)
  dimension :: F_d3 (:,:,:,:,:)
  dimension :: F_p  (:,:,:,:,:)
  dimension :: F_s  (:,:,:,:,:)

  ! local variables ............................................................

  real(RNP), allocatable, save :: v(:,:,:,:,:), w(:,:,:,:,:)
  real(RNP) :: chi_
  logical   :: source_
  integer   :: po, np, ne, nc

  associate( mesh => flow_op % mesh          &
           , bv_u => flow_op % bv_u          &
           , Ms   => flow_op % eop_u % w     &
           , Ds   => flow_op % eop_u % D     &
           , Dd   => flow_op % eop_u % D     )

    ! intialization ............................................................

    po = ubound(Ms,1)
    np = po + 1
    ne = mesh % ne
    nc = problem % nc

    if (present(chi)) then
      chi_ = chi
    else
      chi_ = flow_op % control % chi
    end if

    if (present(source)) then
      source_ = present(F) .and. source
    else
      source_ = present(F)
    end if

    !$omp single
    allocate(w(np,np,np,ne,nc))
    if (present(u_d)) then
      if (problem % HasVariableProperties() .and. .not. present(nu)) then
        allocate(v, mold=w)
      end if
    end if
    !$omp end single

    ! external sources .........................................................

    if (source_ .or. present(F_s)) then
      call problem % GetExternalSources(flow_op%x, t, w)
      if (source_     ) call SetArray(F  , w, multi=.true.)
      if (present(F_s)) call SetArray(F_s, w, multi=.true.)
    else if (present(F)) then
      call SetArray(F, ZERO, multi=.true.)
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
      if (present(F  )) call MergeArrays(ONE , F  , -ONE, w, multi=.true.)
      if (present(F_c)) call MergeArrays(ZERO, F_c, -ONE, w, multi=.true.)

    else if (present(F) .and. present(F_c)) then

      ! use given convection term
      call MergeArrays(ONE, F, ONE, F_c, multi=.true.)
    end if

    ! diffusion ................................................................

    if (present(u_d)) then

      if (present(F_d )) call SetArray(F_d , ZERO, multi = .true.)
      if (present(F_d1)) call SetArray(F_d1, ZERO, multi = .true.)
      if (present(F_d2)) call SetArray(F_d2, ZERO, multi = .true.)
      if (present(F_d3)) call SetArray(F_d3, ZERO, multi = .true.)

      if (present(nu)) then
        call DiffTimeDeriv_Velocity_VI(mesh, Ms, Dd, nu, chi_, u_d, w, F, F_d, &
                                       F_d1, F_d2, F_d3)
        call DiffTimeDeriv_Scalars_VI(mesh, Ms, Dd, nu, u_d, w, F, F_d, F_d1)
      else if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u_d, nu = v)
        call DiffTimeDeriv_Velocity_VI(mesh, Ms, Dd, v, chi_, u_d, w, F, F_d, &
                                       F_d1, F_d2, F_d3)
        call DiffTimeDeriv_Scalars_VI (mesh, Ms, Dd, v, u_d, w, F, F_d, F_d1)
      else
        associate(nu => problem % nu_ref)
          call DiffTimeDeriv_Velocity_CI(mesh, Ms, Dd, nu, chi_, bv_u, u_d, w, &
                                         F, F_d, F_d1, F_d2, F_d3)
          call DiffTimeDeriv_Scalars_CI(mesh, Ms, Dd, nu, u_d, w, F, F_d, F_d1)
        end associate
      end if

    else if (present(F) .and. present(F_d)) then

      ! use given diffusion term
      call MergeArrays(ONE, F, ONE, F_d, multi=.true.)

    end if

    ! pressure .................................................................

    if (present(p)) then
      call WeakGradient(mesh, Ms, Ds, p, w)
      if (present(F)) then
        call MergeArrays(ONE, F(:,:,:,:,:3), -ONE, w(:,:,:,:,:3), multi=.true.)
      end if
      if (present(F_p)) then
        call MergeArrays(ZERO, F_p(:,:,:,:,:3), -ONE, w(:,:,:,:,:3), multi=.true.)
      end if
    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    if (allocated(v)) deallocate(v)
    if (allocated(w)) deallocate(w)
    !$omp end master

  end associate

end subroutine TimeDerivative

!-------------------------------------------------------------------------------
!> Diffusive part of momentum equation with constant viscosity
!>
!> Computes
!>
!>     F(*,1:3) += ν∇²v + ν∇∇·v + χν∇∇·v = F_d(*,1:3)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Velocity_CI(mesh, Ms, Dd, nu, chi, bv_u, v, w, F, F_d, &
    F_d1, F_d2, F_d3)
  class(MeshPartition),    intent(in)    :: mesh            !< mesh partition
  real(RNP),               intent(in)    :: Ms  (:)         !< standard mass matrix
  real(RNP),               intent(in)    :: Dd  (:,:)       !< standard diff matrix
  real(RNP),               intent(in)    :: nu  (:)         !< viscosity
  real(RNP),               intent(in)    :: chi             !< factor of ν∇∇·v
  class(BoundaryVariable), intent(in)    :: bv_u(:)         !< boundary values
  real(RNP),               intent(in)    :: v   (:,:,:,:,:) !< velocity, v(*,1:3)
  real(RNP),               intent(inout) :: w   (:,:,:,:,:) !< workspace
  real(RNP),     optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂v/∂t
  real(RNP),     optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
  real(RNP),     optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇v contribution
  real(RNP),     optional, intent(out)   :: F_d2(:,:,:,:,:) !< ν∇∇·v contribution
  real(RNP),     optional, intent(out)   :: F_d3(:,:,:,:,:) !< χν∇∇·v contribution

  type(BoundaryVariable), allocatable, save :: bv_s(:) ! scalar boundary values
  type(BoundaryVariable), allocatable, save :: bv_v(:) ! vector boundary values

  integer :: ne, np, po
  integer :: b, c

  if (.not. ( present(F)    .or. &
              present(F_d)  .or. &
              present(F_d1) .or. &
              present(F_d2) .or. &
              present(F_d3)     )) return

  ne = mesh%ne
  np = size(Ms)
  po = np - 1

  !$omp single
  allocate(bv_s(mesh % n_boundary))
  allocate(bv_v(mesh % n_boundary))
  do b = 1, mesh % n_boundary
    bv_s(b) = BoundaryVariable(mesh, po, b, nc=1)
    bv_v(b) = BoundaryVariable(mesh, po, b, nc=3)
  end do
  !$omp end single

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    ! div-grad part
    if (present(F) .or. present(F_d) .or. present(F_d1)) then
      do c = 1, 3
        do b = 1, mesh % n_boundary
          bv_s(b) % bc(1)        = bv_u(b) % bc(c)
          bv_s(b) % val(:,:,:,1) = bv_u(b) % val(:,:,:,c)
        end do
       !call TPO_Grad_Eval(np, ne, Dd, mesh%dx, v(:,:,:,:,c), g)
        call WeakGradient(mesh, Ms, Dd, v(:,:,:,:,c), g)
       !call WeakGradient(mesh, Ms, Dd, v(:,:,:,:,c), bv_s, g) !***
        call ScaleArray(g, nu(c), multi=.true.)
       !call WeakDivergence(mesh, Ms, Dd, g, q)
        call TPO_Div_Eval(np, ne, Dd, mesh%dx, g, q)
        if (present(F)) then
          call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
        end if
        if (present(F_d)) then
          call SetArray(F_d(:,:,:,:,c), q)
        end if
        if (present(F_d1)) then
          call SetArray(F_d1(:,:,:,:,c), q)
        end if
      end do
    end if

    ! F_d2 and F_d3
    if (present(F) .or. present(F_d) .or. present(F_d2) .or. present(F_d3)) then
      do b = 1, mesh % n_boundary
        bv_v(b) % bc(1:3)        = bv_u(b) % bc(1:3)
        bv_v(b) % val(:,:,:,1:3) = bv_u(b) % val(:,:,:,1:3)
      end do
     !call TPO_Div_Eval(np, ne, Dd, mesh%dx, v, q)
      call WeakDivergence(mesh, Ms, Dd, v, q)
     !call WeakDivergence(mesh, Ms, Dd, v, bv_v, q) !***
      call ScaleArray(q, sum(nu(1:3))/3)
     !call WeakGradient(mesh, Ms, Dd, q, g)
      call TPO_Grad_Eval(np, ne, Dd, mesh%dx, q, g)
      if (present(F)) then
        call MergeArrays(ONE, F(:,:,:,:,1:3), 1+chi, g, multi=.true.)
      end if
      if (present(F_d)) then
        call MergeArrays(ONE, F_d(:,:,:,:,1:3), 1+chi, g, multi=.true.)
      end if
      if (present(F_d2)) then
        call SetArray(F_d2(:,:,:,:,1:3), g, multi=.true.)
      end if
      if (present(F_d3)) then
        call MergeArrays(ZERO, F_d3(:,:,:,:,1:3), chi, g, multi=.true.)
      end if
    end if

  end associate

  !$omp barrier
  !$omp master
  deallocate(bv_s, bv_v)
  !$omp end master

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

subroutine DiffTimeDeriv_Scalars_CI(mesh, Ms, Dd, nu, u, w, F, F_d, F_d1)
  class(MeshPartition), intent(in)    :: mesh            !< mesh partition
  real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd  (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu  (:)         !< diffusivities
  real(RNP),            intent(in)    :: u   (:,:,:,:,:) !< scalars: u(*,5:nc)
  real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
  real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂u/∂t
  real(RNP),  optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
  real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution

  procedure(TPO_Grad_Proc), pointer :: Gradient

  integer :: nc, ne, np
  integer :: c

  if (.not. (present(F) .or. present(F_d) .or. present(F_d1))) return

  np = size(Ms)
  ne = size(u,4)
  nc = size(u,5)

  call TPO_Grad_Assign(np, Gradient)

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    do c = 5, nc
      call Gradient(np, ne, Dd, mesh%dx, u(:,:,:,:,c), w)
      call ScaleArray(g, nu(c), multi=.true.)
      call WeakDivergence(mesh, Ms, Dd, g, q)
      if (present(F)) then
        call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      end if
      if (present(F_d)) then
        call SetArray(F_d(:,:,:,:,c), q)
      end if
      if (present(F_d1)) then
        call SetArray(F_d1(:,:,:,:,c), q)
      end if
    end do

  end associate

end subroutine DiffTimeDeriv_Scalars_CI

!-------------------------------------------------------------------------------
!> Diffusive part of momentum equation with variable viscosity
!>
!> Computes
!>
!>     F(*,1:3) += ∇·ν[∇v + (∇v)ᵀ] + χ∇ν(∇·v) = F_d(*,1:3)
!>
!> using the weak (projected) form of the outer divergence and gradient
!> operators.

subroutine DiffTimeDeriv_Velocity_VI(mesh, Ms, Dd, nu, chi, v, w, F, F_d, &
    F_d1, F_d2, F_d3)
  class(MeshPartition), intent(in)    :: mesh            !< mesh partition
  real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd  (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu  (:,:,:,:,:) !< viscosity, ν = nu(*,1)
  real(RNP),            intent(in)    :: chi             !< factor of ∇(ν∇·v)
  real(RNP),            intent(in)    :: v   (:,:,:,:,:) !< velocity, v(*,1:3)
  real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
  real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂v/∂t
  real(RNP),  optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
  real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution
  real(RNP),  optional, intent(out)   :: F_d2(:,:,:,:,:) !< ∇·ν(∇v)ᵀ contribution
  real(RNP),  optional, intent(out)   :: F_d3(:,:,:,:,:) !< χ∇ν(∇·v) contribution

  real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:)
  integer :: ne, np
  integer :: c, d, e, i, j, k

  if (.not. ( present(F   ) .or. &
              present(F_d ) .or. &
              present(F_d1) .or. &
              present(F_d2) .or. &
              present(F_d3) ) ) then
    return
  end if

  np = size(Ms)
  ne = mesh%ne

  !$omp single
  allocate(grad_v(np,np,np,ne,3,3))
  !$omp end single

  ! grad_v(c,d) = ∂ v_c / ∂ x_d,  i.e.  grad_v = (∇v)ᵀ
  call TPO_Grad_Eval(np, 3*ne, Dd, mesh%dx, v(:,:,:,:,1:3), grad_v)

  associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

    do c = 1, 3

      ! ∇·ν∇v contribution .....................................................

      if (present(F) .or. present(F_d) .or. present(F_d1)) then

        ! g = ν ∇ v_c
        do d = 1, 3
          !$omp do
          do e = 1, ne
            do k = 1, np
            do j = 1, np
            do i = 1, np
              g(i,j,k,e,d) = nu(i,j,k,e,1) * grad_v(i,j,k,e,c,d)
            end do
            end do
            end do
          end do
        end do

        ! q = ∇·g = ∇·(ν ∇ v_c)
        call WeakDivergence(mesh, Ms, Dd, g, q)

        if (present(F)) then
          call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
        end if
        if (present(F_d )) call SetArray(F_d(:,:,:,:,c), q)
        if (present(F_d1)) call SetArray(F_d1(:,:,:,:,c), q)

      end if

      ! ∇·[ν(∇v)ᵀ] contribution ................................................

      if (present(F) .or. present(F_d) .or. present(F_d2)) then

        ! g = ν (∂ v/ ∂ x_c)
        do d = 1, 3
          !$omp do
          do e = 1, ne
            do k = 1, np
            do j = 1, np
            do i = 1, np
              g(i,j,k,e,d) = nu(i,j,k,e,1) * grad_v(i,j,k,e,d,c)
            end do
            end do
            end do
          end do
        end do

        ! q = ∇·g = ∇·(ν (∂ v/ ∂ x_c))
        call WeakDivergence(mesh, Ms, Dd, g, q)

        if (present(F)) then
          call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
        end if
        if (present(F_d)) then
          call MergeArrays(ONE, F_d(:,:,:,:,c), ONE, q)
        end if
        if (present(F_d2)) then
          call SetArray(F_d2(:,:,:,:,c), q)
        end if

      end if

    end do

    ! div penalty ..............................................................

    if (present(F) .or. present(F_d) .or. present(F_d3)) then

      ! q = ∇·g = χν(∇·v)
      !$omp do
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np
          q(i,j,k,e) = chi * nu(i,j,k,e,1) * ( grad_v(i,j,k,e,1,1) &
                                             + grad_v(i,j,k,e,2,2) &
                                             + grad_v(i,j,k,e,3,3) )
        end do
        end do
        end do
      end do

      ! g = ∇q = χ∇ν(∇·v)
      call WeakGradient(mesh, Ms, Dd, q, g)

      if (present(F)) then
        call MergeArrays(ONE, F(:,:,:,:,1:3), ONE, g, multi = .true.)
      end if
      if (present(F_d)) then
        call MergeArrays(ONE, F_d(:,:,:,:,1:3), ONE, g, multi = .true.)
      end if
      if (present(F_d3)) then
        call SetArray(F_d3(:,:,:,:,1:3), g, multi = .true.)
      end if

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

subroutine DiffTimeDeriv_Scalars_VI(mesh, Ms, Dd, nu, u, w, F, F_d, F_d1)
  class(MeshPartition), intent(in)    :: mesh            !< mesh partition
  real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
  real(RNP),            intent(in)    :: Dd  (:,:)       !< standard diff matrix
  real(RNP),            intent(in)    :: nu  (:,:,:,:,:) !< diffusivities
  real(RNP),            intent(in)    :: u   (:,:,:,:,:) !< scalars: u(*,5:nc)
  real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
  real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂u/∂t
  real(RNP),  optional, intent(inout) :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
  real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution

  procedure(TPO_Grad_Proc), pointer :: Gradient

  integer :: nc, ne, np
  integer :: c, d, e, i, j, k

  if (.not. (present(F) .or. present(F_d) .or. present(F_d1))) return

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

      if (present(F)) then
        call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
      end if
      if (present(F_d )) call SetArray(F_d (:,:,:,:,c), q)
      if (present(F_d1)) call SetArray(F_d1(:,:,:,:,c), q)

    end do

  end associate

end subroutine DiffTimeDeriv_Scalars_VI

!===============================================================================

end module CART__ISP_Flow__Time_Derivative
