!> summary:  ISP flow: computation of time derivative to given solution
!> author:   Joerg Stiller
!> date:     2018/05/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: computation of time derivative to given solution
!===============================================================================

module CART__ISP_Flow__Time_Derivative

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Array_Assignments
  use TPO_sDDD
  use ISP_Flow_Problem
  use CART__TPO_Div
  use CART__TPO_Grad
  use CART__TPO_RotRot
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
  integer :: c

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
      call AssignArray(F_s, F, multi=.true.)
    end if

    ! convection ...............................................................

    if (problem%stokes) then

      ! skip convection term
      if (present(F_c)) then
        call AssignScalar(F_c, ZERO, multi=.true.)
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

    if (present(u_d)) then ! so far assuming constant nu

      associate(v => u_d(:,:,:,:,1:3), q => w(:,:,:,:,4), nu => problem%nu_ref)

        ! div-grad part
        do c = 1, nc
          if (c == 4) then
            call AssignScalar(q, ZERO)
          else
            call TPO_Grad_Eval(np, ne, Dd, mesh%dx, u_d(:,:,:,:,c), w)
            call ScaleArray(w, nu(c), multi=.true.)
            call WeakDivergence(mesh, Ms, Dd, w, q)
            call MergeArrays(ONE, F(:,:,:,:,c), ONE, q)
          end if
          if (present(F_d)) then
            call AssignArray(F_d(:,:,:,:,c), q)
          end if
        end do

        ! divergence penalty
        call TPO_Div_Eval(np, ne, Ds, mesh%dx, v, q)
        call ScaleArray(q, -nu(1))
        call WeakGradient(mesh, Ms, Ds, q, w)
        call MergeArrays(ONE, F, ONE, w, multi=.true.)
        if (present(F_d)) then
          call MergeArrays(ONE, F_d, ONE, w, multi=.true.)
        end if

      end associate

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

!===============================================================================

end module CART__ISP_Flow__Time_Derivative