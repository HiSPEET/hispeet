!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  ISP flow: projection step
!> author:   Joerg Stiller
!> date:     2019/05/08
!>
!> @todo
!>   * put DivergencePenalty into separate module
!>   * use combined grad-div TPO instead of subsequent calls to div and grad
!>   * test using velocity BV when computing the traces
!>   * check element-based assembly of jump contributions
!> @endtodo
!===============================================================================

module CART__ISP_Flow__Projection

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, TWO, THIRD
  use XMPI

  use Array_Assignments
  use Array_Reductions
  use Standard_Element_Operators__1D
  use TPO__Diagonal__3D
  use TPO__Div__3D_R
  use TPO__Grad__3D_R

  use ISP_Flow_Problem

  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Normal_Trace_Operator
  use CART__DG_Weak_Gradient
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: ProjectionStep
  public :: DivergencePenalty
  public :: DivergencePenaltyFactor

contains

  !-----------------------------------------------------------------------------
  !>

  subroutine ProjectionStep(problem, flow_op, dt, p, u, w)
    class(FlowProblem),   intent(in)    :: problem      !< flow problem
    class(FlowOperators), intent(in)    :: flow_op      !< flow operators
    real(RNP),            intent(in)    :: dt           !< time step size
    real(RNP),            intent(in)    :: p(:,:,:,:)   !< pressure
    real(RNP),            intent(inout) :: u(:,:,:,:,:) !< current solution
    real(RNP),            intent(inout) :: w(:,:,:,:,:) !< workspace

    real(RNP) :: gamma

    associate( mesh => flow_op % mesh     &
             , eop  => flow_op % eop_u    &
             , v    => u(:,:,:,:,1:3)     &
             , f    => w(:,:,:,:,1:3)     )

      ! pressure correction ....................................................

      ! v = v - ∆t ∇p
      call WeakGradient(mesh, eop%w, eop%D, p, f)    ! f = ∇p
      call MergeArrays(ONE, v, -dt, f, multi=.true.) ! v = v - ∆t f

      ! divergence penalization ................................................

      gamma = DivergencePenaltyFactor( penalty = flow_op % control % div_penalty &
                                     , coeff   = flow_op % control % div_coeff   &
                                     , nu      = problem % nu_ref(1)             &
                                     , dt      = dt                              &
                                     , dx      = mesh % dx                       &
                                     , po      = eop % po                        )

      if (gamma > 0) then
        call SetArray(f, v, multi = .true.)
        call ProjectionSolver( mesh, eop, dt, gamma           &
                             , flow_op % control % div_r_red  &
                             , flow_op % control % div_i_max  &
                             , f, v                           )
      end if

    end associate

  end subroutine ProjectionStep

!-----------------------------------------------------------------------------
!>

subroutine ProjectionSolver(mesh, sop, dt, gamma, r_red, i_max, f, v)
  class(MeshPartition),               intent(in) :: mesh !< mesh partition
  class(StandardElementOperators_1D), intent(in) :: sop  !< std SE operators
  real(RNP), intent(in)    :: dt           !< time step size
  real(RNP), intent(in)    :: gamma        !< penalty factor
  real(RNP), intent(in)    :: r_red        !< min res reduction
  integer,   intent(in)    :: i_max        !< max num iterations
  real(RNP), intent(in)    :: f(:,:,:,:,:) !< RHS
  real(RNP), intent(inout) :: v(:,:,:,:,:) !< velocity

  ! local variables ..........................................................

  real(RNP), dimension(:,:,:,:,:), allocatable, save :: r, p, q, z
  logical  , save :: converged

  real(RNP), allocatable :: Ms_inv(:)
  real(RNP) :: alpha, beta, rz, rz_old, rz_max
  real(RNP) :: g, g_inv
  integer   :: np, ne
  integer   :: i, k

  ! initialization .............................................................

  associate(Ms => sop%w, Ds => sop%D)

    np = size(v,1)
    ne = mesh % ne

    ! work space
    !$omp single
    allocate(r, mold = v)
    allocate(p, mold = v)
    allocate(q, mold = v)
    allocate(z, mold = v)
    !$omp end single

    allocate(Ms_inv, source = 1/Ms)

    g = product(mesh%dx) / 8
    g_inv = 1 / g

    ! initial residual .........................................................

    ! z₀ = f - (v + Δt DivPenalty(v))
    call DivergencePenalty(mesh, Ms, Ds, gamma, v, z)  ! z₀ = DivPenalty(v)
    call MergeArrays(  dt, z, ONE, v, multi = .true.)  ! z₀ = Δt z₀ + v
    call MergeArrays(-ONE, z, ONE, f, multi = .true.)  ! z₀ = f  - z₀

    ! r₀ = M z₀
    do k = 1, 3
      call TPO_Diagonal(g, Ms, z(:,:,:,:,k), r(:,:,:,:,k))
    end do

    ! residual metrics
    rz = ScalarProduct(r, z, mesh%comm)
    rz_max = rz * r_red**2 / g

    ! initial search vector
    call SetArray(p, z, multi = .true.)

    do i = 1, i_max

      rz_old = rz

      ! q = M (p + Δt DivPenalty(p))
      call DivergencePenalty(mesh, Ms, Ds, gamma, p, z) ! z = DivPenalty(p)
      call MergeArrays(dt, z, ONE, p, multi = .true.)   ! z = Δt z + p
      do k = 1, 3
        call TPO_Diagonal(g, Ms, z(:,:,:,:,k), q(:,:,:,:,k)) ! q = M z
      end do

      alpha = rz / ScalarProduct(p, q, mesh%comm)

      ! correction
      call MergeArrays(ONE, v,  alpha, p, multi = .true.) ! v = v + alpha p
      call MergeArrays(ONE, r, -alpha, q, multi = .true.) ! r = r - alpha q
      do k = 1, 3
        call TPO_Diagonal(g_inv, Ms_inv, r(:,:,:,:,k), z(:,:,:,:,k)) ! z = M⁻¹ r
      end do

      rz = ScalarProduct(r, z, mesh%comm)

      !$omp master
      converged = rz <= rz_max
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier

      if (converged .or. i == i_max) exit

      ! new search vector
      beta = rz / rz_old
      call MergeArrays(beta, p,  ONE, z, multi = .true.)  ! p = z + beta p

    end do

    ! finalization .............................................................

    !$omp barrier
    !$omp master
    deallocate(r, p, q, z)
    !$omp end master

  end associate

end subroutine ProjectionSolver

!-------------------------------------------------------------------------------
!> Divergence penalty factor

function DivergencePenaltyFactor(penalty, coeff, nu, dt, dx, po) result(gamma)
  integer,   intent(in) :: penalty !< type of penalty
  real(RNP), intent(in) :: coeff   !< scaling coefficient
  real(RNP), intent(in) :: nu      !< viscosity
  real(RNP), intent(in) :: dt      !< time step size
  real(RNP), intent(in) :: dx(3)   !< element extensions
  integer,   intent(in) :: po      !< polynomial order
  real(RNP)             :: gamma

  real(RNP) :: h

  select case(penalty)
  case(1)
    gamma = coeff * nu
  case(2)
    h = product(dx) ** THIRD
    gamma = coeff * (h/po)**2 / dt
  case default
    gamma = 0
  end select

end function DivergencePenaltyFactor

!-------------------------------------------------------------------------------
!> Divergence penalty scaled with the inverse mass matrix
!>
!> Computes
!>
!>     w = M⁻¹γ (∫(∇·φ)(∇·v)dΩ + 1/hᵉ ∫[φ][v]dΓ)
!>
!> where φ represents the vector-valued test functions constructed from single
!> base functions, v the given velocity, [φ] and [v] their respective jumps
!> across the faces and hᵉ the average length of the element supporting φ

subroutine DivergencePenalty(mesh, Ms, Ds, gamma, v, w)
  class(MeshPartition), intent(in)  :: mesh            !< mesh partition
  real(RNP),            intent(in)  :: Ms(0:)          !< std mass matrix
  real(RNP),            intent(in)  :: Ds(0:,0:)       !< std diff matrix
  real(RNP),            intent(in)  :: gamma           !< penalty factor
  real(RNP),            intent(in)  :: v(0:,0:,0:,:,:) !< velocity field
  real(RNP),            intent(out) :: w(0:,0:,0:,:,:) !< div-penalty term

  ! internal variables .........................................................

  type(NormalTraceOperator), asynchronous, allocatable, save :: trace_op
  real(RNP), allocatable, save :: tr_vn(:,:,:,:)
  real(RNP), allocatable, save :: div_v(:,:,:,:)

  real(RNP), allocatable :: Dt(:,:), wf(:,:)
  real(RNP) :: g(3)

  integer :: po, ne, np
  integer :: tag = 1000
  integer :: e, i, j, f, f1, f2

  ! initialization .............................................................

  po = ubound(Ms,1)
  ne = mesh % ne
  np = po + 1

  !$omp single
  allocate(trace_op)
  allocate(tr_vn(0:po, 0:po, 2, mesh%nf))
  allocate(div_v(0:po, 0:po, 0:po, ne))
  !$omp end single

  ! extract and start transferring traces ......................................

  ! NB: it may be beneficial to include the velocity BC here

  call trace_op % GetTrace_Start(mesh, v, tr_vn, tag)

  ! element-interior contributions  ............................................

  ! NB: this should be replaced by a single grad-div TPO

  ! transposed element gradient operator
  allocate(Dt(0:po,0:po))
  do j = 0, po
  do i = 0, po
    Dt(i,j) = Ms(j) * Ds(j,i) / Ms(i)
  end do
  end do

  ! divergence
  call TPO_Div_R(Ms, Ds, mesh%dx, u=v, v=div_v)

  ! transposed gradient with zero boundary flux, including penalty factor
  call TPO_Grad_R(Ms, Dt, mesh%dx, div_v, v=w)
  call ScaleArray(w(:,:,:,:,1:3), gamma)

  ! element-boundary contributions  ............................................

  allocate(wf(0:po,0:po))

  ! complete transfer
  call trace_op % GetTrace_Finish(mesh, tr_vn)

  ! metric including penalty
  g = 2 / (Ms(0) * mesh%dx) * gamma

  associate(dx => mesh%dx, face => mesh%face)

    ! x1-flux contribution .....................................................

    f1 = 1                 ! first face
    f2 = mesh % nf1        ! last

    !$omp do private(wf)
    do f = f1, f2

      wf = g(1) * (tr_vn(:,:,1,f) + tr_vn(:,:,2,f))

      e = face(f) % element(1)
      if (0 < e .and. e <= mesh % ne) then
        w(po,:,:,e,1) = w(po,:,:,e,1) + wf
      end if

      e = face(f) % element(2)
      if (0 < e .and. e <= mesh % ne) then
        w( 0,:,:,e,1) = w( 0,:,:,e,1) - wf
      end if

    end do

    ! x2-flux contribution .....................................................

    f1 = f2 + 1            ! first face
    f2 = f2 + mesh % nf2   ! last

    !$omp do private(wf)
    do f = f1, f2

      wf = g(2) * (tr_vn(:,:,1,f) + tr_vn(:,:,2,f))

      e = face(f) % element(1)
      if (0 < e .and. e <= mesh % ne) then
        w(:,po,:,e,2) = w(:,po,:,e,2) + wf
      end if

      e = face(f) % element(2)
      if (0 < e .and. e <= mesh % ne) then
        w(:, 0,:,e,2) = w(:, 0,:,e,2) - wf
      end if

    end do

    ! x3-flux contribution .....................................................

    f1 = f2 + 1            ! first face
    f2 = f2 + mesh % nf3   ! last

    !$omp do private(wf)
    do f = f1, f2

      wf = g(3) * (tr_vn(:,:,1,f) + tr_vn(:,:,2,f))

      e = face(f) % element(1)
      if (0 < e .and. e <= mesh % ne) then
        w(:,:,po,e,3) = w(:,:,po,e,3) + wf
      end if

      e = face(f) % element(2)
      if (0 < e .and. e <= mesh % ne) then
        w(:,:,0,e,3) = w(:,:,0,e,3) - wf
      end if

    end do

  end associate

  ! finalization ...............................................................

  !$omp barrier
  !$omp master
  deallocate(tr_vn)
  deallocate(trace_op)
  deallocate(div_v)
  !$omp end master

end subroutine DivergencePenalty

!===============================================================================

end module CART__ISP_Flow__Projection
