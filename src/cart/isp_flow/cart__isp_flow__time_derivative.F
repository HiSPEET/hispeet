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
!>      F_c = -∇·(v u)
!>
!> with `v = u(*,1:3)`,
!>
!>      F_d = F_d1 + F_d2 + F_d3
!>
!>      F_d1(*, : ) = ∇·(ν (∇v) )
!>      F_d2(*,1:3) = ∇·(ν (∇v)ᵀ)
!>      F_d3(*,1:3) = χ∇(ν (∇·v))
!>
!> with `v = u(*,1:3)`,
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
  use TPO__Div__3D_R
  use TPO__Grad__3D_R
  use ISP_Flow_Problem
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

  !-----------------------------------------------------------------------------
  !> Computation of time derivative to given solution
  !>
  !> @note
  !>   * if absent, the bulk diffusion parameter χ ist taken from `flow_op`
  !> @endnote

  subroutine TimeDerivative( problem, flow_op, t, u, p, nu, chi,             &
                             F, F_c, F_d, F_d1, F_d2, F_d3, F_p, F_s, source )

    ! arguments ................................................................

    class(FlowProblem),   intent(in)  :: problem !< flow problem
    class(FlowOperators), intent(in)  :: flow_op !< flow operators
    real(RNP),            intent(in)  :: t       !< time
    real(RNP),            intent(in)  :: u       !< flow variables
    real(RNP),  optional, intent(in)  :: p       !< pressure
    real(RNP),  optional, intent(in)  :: nu      !< variable diffusivity
    real(RNP),  optional, intent(in)  :: chi     !< bulk diffusion parameter χ
    real(RNP),  optional, intent(out) :: F       !< ∂u/∂t
    real(RNP),  optional, intent(out) :: F_c     !< convection part
    real(RNP),  optional, intent(out) :: F_d     !< diffusion part, complete
    real(RNP),  optional, intent(out) :: F_d1    !< diffusion, only ∇·ν∇u
    real(RNP),  optional, intent(out) :: F_d2    !< diffusion, only ∇·ν(∇v)ᵀ
    real(RNP),  optional, intent(out) :: F_d3    !< diffusion, only χ∇ν(∇·v)
    real(RNP),  optional, intent(out) :: F_p     !< pressure part
    real(RNP),  optional, intent(out) :: F_s     !< source part
    logical,    optional, intent(in)  :: source  !< include sources in F [T]

    dimension :: u    (:,:,:,:,:)
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

    ! local variables ..........................................................

    real(RNP), allocatable, save :: v(:,:,:,:,:), w(:,:,:,:,:)
    real(RNP) :: chi_
    logical   :: source_
    integer   :: po, np, ne, nc

    ! with SVV
    real(RNP), allocatable :: Ds_svv(:,:) ! standard SVV diff matrix

    associate( mesh => flow_op % mesh          &
             , bv_u => flow_op % bv_u          &
             , Ms   => flow_op % eop_u % w     &
             , Ds   => flow_op % eop_u % D     )

      ! intialization ..........................................................

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
      allocate(v(np,np,np,ne,nc))
      allocate(w(np,np,np,ne,nc))
      !$omp end single

      ! external sources .......................................................

      if (source_ .or. present(F_s)) then
        call problem % GetExternalSources(flow_op%x, t, w)
        if (source_     ) call SetArray(F  , w, multi=.true.)
        if (present(F_s)) call SetArray(F_s, w, multi=.true.)
      else if (present(F)) then
        call SetArray(F, ZERO, multi=.true.)
      end if

      ! convection .............................................................

      if (problem%stokes) then

        ! skip convection term
        if (present(F_c)) then
          call SetArray(F_c, ZERO, multi=.true.)
        end if

      else if (present(F) .or. present(F_c)) then

        ! compute convection term
        call WeakConvectiveFlux(flow_op, u, div_F=w)
        if (present(F  )) call MergeArrays(ONE , F  , -ONE, w, multi=.true.)
        if (present(F_c)) call MergeArrays(ZERO, F_c, -ONE, w, multi=.true.)

      end if

      ! diffusion ..............................................................

      if (present(F_d )) call SetArray(F_d (:,:,:,:,4:nc), ZERO, multi = .true.)
      if (present(F_d1)) call SetArray(F_d1(:,:,:,:,4   ), ZERO                )
      if (present(F_d2)) call SetArray(F_d2(:,:,:,:,4:nc), ZERO, multi = .true.)
      if (present(F_d3)) call SetArray(F_d3(:,:,:,:,4:nc), ZERO, multi = .true.)

      if (present(nu)) then

        ! variable diffusivity nu given
        call DiffTimeDeriv_Velocity_VI( mesh, Ms, Ds, nu, chi_, u, w, &
                                        F, F_d, F_d1, F_d2, F_d3      )
        call DiffTimeDeriv_Scalars_VI ( mesh, Ms, Ds, nu, u, w, F, F_d, F_d1 )

      else if (problem % HasVariableProperties()) then

        ! variable diffusivity nu absent
        call problem % GetDiffusivity(flow_op%x, t, u, nu = v)
        call DiffTimeDeriv_Velocity_VI( mesh, Ms, Ds, v, chi_, u, w, &
                                        F, F_d, F_d1, F_d2, F_d3     )
        call DiffTimeDeriv_Scalars_VI ( mesh, Ms, Ds, v, u, w, F, F_d, F_d1 )

      else if (flow_op % eop_u % Has_SVV()) then

        ! constant isotropic viscosity and SVV
        associate( nu     => problem % nu_ref,    &
                   nu_svv => problem % nu_svv_ref )

          allocate(Ds_svv(0:po, 0:po))
          call flow_op % eop_u % Get_SVV_StandardDiffMatrix(Ds_svv(:,:))
          call DiffTimeDeriv_Velocity_CI_svv( mesh, Ms, Ds, Ds_svv     &
                                            , nu, nu_svv, chi_, u, w   &
                                            , F, F_d, F_d1, F_d2, F_d3 )
         call DiffTimeDeriv_Scalars_CI_svv ( mesh, Ms, Ds, Ds_svv     &
                                            , nu, nu_svv, u, w         &
                                            , F, F_d1                  )
        end associate

      else

        ! constant isotropic viscosity and no SVV
        associate(nu => problem % nu_ref)
          call DiffTimeDeriv_Velocity_CI( mesh, Ms, Ds, nu, chi_, u, w,  &
                                          F, F_d, F_d1, F_d2, F_d3)
          call DiffTimeDeriv_Scalars_CI ( mesh, Ms, Ds, nu, u, w, F, F_d, F_d1 )
        end associate

      end if

      ! pressure ...............................................................

      if (present(p)) then
        call WeakGradient(mesh, Ms, Ds, p, w)
        if (present(F)) then
          call MergeArrays(ONE, F(:,:,:,:,:3), -ONE, w(:,:,:,:,:3), multi=.true.)
        end if
        if (present(F_p)) then
          call MergeArrays(ZERO, F_p(:,:,:,:,:3), -ONE, w(:,:,:,:,:3), multi=.true.)
        end if
      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp single
      if (allocated(v)) deallocate(v)
      if (allocated(w)) deallocate(w)
      !$omp end single

    end associate

  end subroutine TimeDerivative

  !-----------------------------------------------------------------------------
  !> Diffusive part of momentum equation with constant viscosity
  !>
  !> Computes
  !>
  !>     F(*,1:3) += ν∇²v + ν∇∇·v + χν∇∇·v = F_d(*,1:3)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Velocity_CI( mesh, Ms, Ds, nu, chi, v, w &
                                      , F, F_d, F_d1, F_d2, F_d3    )
    class(MeshPartition),    intent(in)    :: mesh            !< mesh partition
    real(RNP),               intent(in)    :: Ms  (:)         !< standard mass matrix
    real(RNP),               intent(in)    :: Ds  (:,:)       !< standard diff matrix
    real(RNP),               intent(in)    :: nu  (:)         !< viscosity
    real(RNP),               intent(in)    :: chi             !< factor of ν∇∇·v
    real(RNP),               intent(in)    :: v   (:,:,:,:,:) !< velocity, v(*,1:3)
    real(RNP),               intent(inout) :: w   (:,:,:,:,:) !< workspace
    real(RNP),     optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂v/∂t
    real(RNP),     optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
    real(RNP),     optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇v contribution
    real(RNP),     optional, intent(out)   :: F_d2(:,:,:,:,:) !< ν∇∇·v contribution
    real(RNP),     optional, intent(out)   :: F_d3(:,:,:,:,:) !< χν∇∇·v contribution

    integer :: ne, np, po
    integer :: c

    if (.not. ( present(F)    .or. &
                present(F_d)  .or. &
                present(F_d1) .or. &
                present(F_d2) .or. &
                present(F_d3)     )) return

    ne = mesh%ne
    np = size(Ms)
    po = np - 1

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      ! div-grad part
      if (present(F) .or. present(F_d) .or. present(F_d1)) then
        do c = 1, 3
          call WeakGradient(mesh, Ms, Ds, v(:,:,:,:,c), g)
          call TPO_Div_R(Ms, nu(c) * Ds, mesh%dx, u=g, v=q)
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
        call WeakDivergence(mesh, Ms, Ds, v, q)
        call TPO_Grad_R(Ms, sum(nu(1:3))/3 * Ds, mesh%dx, q, v=g)
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

  end subroutine DiffTimeDeriv_Velocity_CI

  !-----------------------------------------------------------------------------
  !> Diffusive part of scalar transport with constant diffusivity
  !>
  !> Computes
  !>
  !>     F(*,5:nc) += ν∇·∇u = F_d(*,5:nc)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Scalars_CI(mesh, Ms, Ds, nu, u, w, F, F_d, F_d1)
    class(MeshPartition), intent(in)    :: mesh            !< mesh partition
    real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
    real(RNP),            intent(in)    :: Ds  (:,:)       !< standard diff matrix
    real(RNP),            intent(in)    :: nu  (:)         !< diffusivities
    real(RNP),            intent(in)    :: u   (:,:,:,:,:) !< scalars: u(*,5:nc)
    real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
    real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂u/∂t
    real(RNP),  optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
    real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution

    integer :: nc, ne, np
    integer :: c

    if (.not. (present(F) .or. present(F_d) .or. present(F_d1))) return

    np = size(Ms)
    ne = size(u,4)
    nc = size(u,5)

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      do c = 5, nc
        call TPO_Grad_R(Ms, nu(c) * Ds, mesh%dx, u(:,:,:,:,c), v=w)
        call WeakDivergence(mesh, Ms, Ds, g, q)
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
  !> Diffusive part of momentum equation with a combination of constant viscosity
  !> and constant spectral diffusivity
  !> Computes
  !>
  !>     F(*,1:3) += ν∇·[∇v + (∇v)ᵀ] + χν∇∇·v + ∇·(νˢQ∇v) = F_d(*,1:3)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Velocity_CI_svv( mesh, Ms, Ds, Ds_svv     &
                                          , nu, nu_svv, chi, v, w    &
                                          , F, F_d, F_d1, F_d2, F_d3 )

    class(MeshPartition), intent(in)    :: mesh            !< mesh partition
    real(RNP),            intent(in)    :: Ms(:)           !< standard mass matrix
    real(RNP),            intent(in)    :: Ds(:,:)         !< standard diff matrix
    real(RNP),            intent(in)    :: Ds_svv(:,:)     !< std SVV diff matrix
    real(RNP),            intent(in)    :: nu(:)           !< viscosity
    real(RNP),            intent(in)    :: nu_svv(:)       !< SVV
    real(RNP),            intent(in)    :: chi             !< factor of ν∇∇·v
    real(RNP),            intent(in)    :: v   (:,:,:,:,:) !< velocity, v(*,1:3)
    real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
    real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂v/∂t
    real(RNP),  optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
    real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·(ν∇v+νˢQ∇v) contrib.
    real(RNP),  optional, intent(out)   :: F_d2(:,:,:,:,:) !< ν∇∇·v contribution
    real(RNP),  optional, intent(out)   :: F_d3(:,:,:,:,:) !< χν∇∇·v contribution

    real(RNP), allocatable :: Bs(:,:) ! standard viscous flux matrix
    integer :: ne, np
    integer :: c

    if (.not. ( present(F)    .or. &
                present(F_d)  .or. &
                present(F_d1) .or. &
                present(F_d2) .or. &
                present(F_d3)     )) return

    np = size(Ms)
    ne = mesh%ne

    allocate(Bs, mold = Ds)

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      ! div-grad part
      if (present(F) .or. present(F_d) .or. present(F_d1)) then
        do c = 1, 3
          Bs = nu(c) * Ds + nu_svv(c) * Ds_svv
          ! computes the gradient scaled with the viscosity ν∇v + νˢQ∇v
          call TPO_Grad_R(Ms, Bs, mesh%dx, v(:,:,:,:,c), v = g)
          call WeakDivergence(mesh, Ms, Ds, g, q)
          if (present(F)) then
            call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
          end if
          if (present(F_d)) then
            call SetArray(F_d (:,:,:,:,c), q)
          end if
          if (present(F_d1)) then
            call SetArray(F_d1(:,:,:,:,c), q)
          end if
        end do
      end if

      ! F_d2 and F_d3
      if (present(F) .or. present(F_d) .or. present(F_d2) .or. present(F_d3)) then
        call WeakDivergence(mesh, Ms, Ds, v, q)
        call TPO_Grad_R(Ms, sum(nu(1:3))/3 * Ds, mesh%dx, q, v = g)
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

  end subroutine DiffTimeDeriv_Velocity_CI_svv

  !-----------------------------------------------------------------------------
  !> Diffusive part of scalar transport with constant diffusivity
  !>
  !> Computes
  !>
  !>     F(*,5:nc) += ∇·(ν∇u + νˢQ∇u) = F_d(*,5:nc)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Scalars_CI_svv( mesh, Ms, Ds, Ds_svv, nu, nu_svv &
                                         , u, w, F, F_d, F_d1               )
    class(MeshPartition), intent(in)    :: mesh            !< mesh partition
    real(RNP),            intent(in)    :: Ms(:)           !< standard mass matrix
    real(RNP),            intent(in)    :: Ds(:,:)         !< standard diff matrix
    real(RNP),            intent(in)    :: Ds_svv(:,:)     !< std SVV diff matrix
    real(RNP),            intent(in)    :: nu(:)           !< viscosity
    real(RNP),            intent(in)    :: nu_svv(:)       !< SVV
    real(RNP),            intent(in)    :: u   (:,:,:,:,:) !< scalars: u(*,5:nc)
    real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
    real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂u/∂t
    real(RNP),  optional, intent(out)   :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
    real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution

    real(RNP), allocatable :: Bs(:,:) ! standard viscous flux matrix
    integer :: nc, ne, np
    integer :: c

    if (.not. (present(F) .or. present(F_d) .or. present(F_d1))) return

    np = size(Ms)
    ne = size(u,4)
    nc = size(u,5)

    allocate(Bs, mold = Ds)

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      do c = 5, nc
        Bs = nu(c) * Ds + nu_svv(c) * Ds_svv
        ! computes the gradient scaled with the viscosity ν∇v + νˢQ∇v
        call TPO_Grad_R(Ms, Bs(:,:), mesh%dx, u(:,:,:,:,c), v = w)
        call WeakDivergence(mesh, Ms, Ds, g, q)
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

  end subroutine DiffTimeDeriv_Scalars_CI_svv

  !-----------------------------------------------------------------------------
  !> Diffusive part of momentum equation with variable viscosity
  !>
  !> Computes
  !>
  !>     F(*,1:3) += ∇·ν[∇v + (∇v)ᵀ] + χ∇ν(∇·v) = F_d(*,1:3)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Velocity_VI(mesh, Ms, Ds, nu, chi, v, w, F, F_d, &
      F_d1, F_d2, F_d3)
    class(MeshPartition), intent(in)    :: mesh            !< mesh partition
    real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
    real(RNP),            intent(in)    :: Ds  (:,:)       !< standard diff matrix
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

    ! grad_v(d,c) = ∂ v_c / ∂ x_d,  i.e.  grad_v = ∇v
    do c = 1, 3
      call TPO_Grad_R(Ms, Ds, mesh%dx, v(:,:,:,:,c), v = grad_v(:,:,:,:,:,c))
    end do

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      do c = 1, 3

        ! ∇·ν∇v contribution ...................................................

        if (present(F) .or. present(F_d) .or. present(F_d1)) then

          ! g = ν ∇ v_c
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

          ! q = ∇·g = ∇·(ν ∇ v_c)
          call WeakDivergence(mesh, Ms, Ds, g, q)

          if (present(F)) then
            call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
          end if
          if (present(F_d )) call SetArray(F_d(:,:,:,:,c), q)
          if (present(F_d1)) call SetArray(F_d1(:,:,:,:,c), q)

        end if

        ! ∇·[ν(∇v)ᵀ] contribution ..............................................

        if (present(F) .or. present(F_d) .or. present(F_d2)) then

          ! g = ν (∂ v/ ∂ x_c)
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

          ! q = ∇·g = ∇·(ν (∂ v/ ∂ x_c))
          call WeakDivergence(mesh, Ms, Ds, g, q)

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

      ! div penalty ............................................................

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
        call WeakGradient(mesh, Ms, Ds, q, g)

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

  !-----------------------------------------------------------------------------
  !> Diffusive part of scalar transport with variable diffusivity
  !>
  !> Computes
  !>
  !>     F(*,5:nc) += ∇·ν∇u = F_d(*,5:nc)
  !>
  !> using the weak (projected) form of the outer divergence and gradient
  !> operators.

  subroutine DiffTimeDeriv_Scalars_VI(mesh, Ms, Ds, nu, u, w, F, F_d, F_d1)
    class(MeshPartition), intent(in)    :: mesh            !< mesh partition
    real(RNP),            intent(in)    :: Ms  (:)         !< standard mass matrix
    real(RNP),            intent(in)    :: Ds  (:,:)       !< standard diff matrix
    real(RNP),            intent(in)    :: nu  (:,:,:,:,:) !< diffusivities
    real(RNP),            intent(in)    :: u   (:,:,:,:,:) !< scalars: u(*,5:nc)
    real(RNP),            intent(inout) :: w   (:,:,:,:,:) !< workspace
    real(RNP),  optional, intent(inout) :: F   (:,:,:,:,:) !< time derivative ∂u/∂t
    real(RNP),  optional, intent(inout) :: F_d (:,:,:,:,:) !< ∂u/∂t diffusive part
    real(RNP),  optional, intent(out)   :: F_d1(:,:,:,:,:) !< ∇·ν∇u contribution

    integer :: nc, ne, np
    integer :: c, d, e, i, j, k

    if (.not. (present(F) .or. present(F_d) .or. present(F_d1))) return

    np = size(Ms)
    ne = size(u,4)
    nc = size(u,5)

    associate(g => w(:,:,:,:,1:3), q => w(:,:,:,:,4))

      do c = 5, nc

        call TPO_Grad_R(Ms, Ds, mesh%dx, u(:,:,:,:,c), v=g)

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

        call WeakDivergence(mesh, Ms, Ds, g, q)

        if (present(F)) then
          call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
        end if
        if (present(F_d )) call SetArray(F_d (:,:,:,:,c), q)
        if (present(F_d1)) call SetArray(F_d1(:,:,:,:,c), q)

      end do

    end associate

  end subroutine DiffTimeDeriv_Scalars_VI

  !=============================================================================

end module CART__ISP_Flow__Time_Derivative
