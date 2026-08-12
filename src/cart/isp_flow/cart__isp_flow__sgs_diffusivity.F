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

!> summary:  SGS diffusivity for Cartesian incompressible flow solver
!> author:   Yang Liu, Martina Grotteschi, Joerg Stiller
!> date:     2020/11/06
!>
!> @note
!>   Placement of dynamic procedure in element-based routine problematic
!>     - inefficient memory usage
!>     - allocation and deallocation of identical workspace for each element
!>     - major revision required for parallelization with OpenMP
!>
!===============================================================================

module CART__ISP_Flow__SGS_Diffusivity

  use Kind_Parameters, only: RNP
  use Constants      , only: ZERO, ONE, HALF, THIRD
  use Array_Assignments

  use ISP_Flow_SGS__Sigma_Model

  use TPO__Grad__3D_R         ! element-based computation of grad(u)
  use CART__DG_Weak_Gradient  ! L2 projection of grad(u) [overkill?]
  use CART__Mesh_Partition
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: Add_SGS_Diffusivity, GetDynamicCoefficient

contains

  !-----------------------------------------------------------------------------
  !> Subgrid-scale turbulence diffusivity
  !>
  !> Evaluates SGS diffusivity and adds it to the given physical diffusivity

  subroutine Add_SGS_Diffusivity(flow_op, u, nu)
    class(FlowOperators),  intent(in)    :: flow_op       !< flow operators
    real(RNP), contiguous, intent(in)    :: u(:,:,:,:,:)  !< flow variables
    real(RNP), contiguous, intent(inout) :: nu(:,:,:,:,:) !< diffusivity

    real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:) ! velocity gradient
    real(RNP), allocatable, save :: c_sigma(:)          ! model constant

    integer :: np, ne
    integer :: i

    !$omp master

    np = size(u,1) ! number of element points per direction
    ne = size(u,4) ! number of elements

    associate( mesh      => flow_op % mesh  & ! access to mesh
             , eop       => flow_op % eop_u & ! access to element operators
             , sgs_coeff => flow_op % control % sgs_coeff ) ! SGS model constant

      ! velocity gradient ......................................................

      allocate(grad_v(np,np,np,ne,1:3,1:3))
      allocate(c_sigma(ne))

      ! compute velocity gradient, grad_v(:,:,:,:,1:3,i) = ∇vᵢ
      do i = 1, 3
        call TPO_Grad_R(eop%w, eop%D, mesh%dx, u(:,:,:,:,i), v=grad_v(:,:,:,:,:,i))
      end do

      ! evaluation of model constant ...........................................

      if (flow_op % control % sgs_dyn) then
        ! evaluate the dynamic model constant
        call GetDynamicCoefficient(flow_op, u, grad_v, c_sigma)
      else
        ! constant model constant
        do i = 1, ne
          c_sigma(i) = sgs_coeff
        end do
      end if

      ! add nu evaluated using sigma model
      call EvalSigmaModel(c_sigma, mesh%dx, grad_v, nu)

      deallocate(grad_v)
      deallocate(c_sigma)

    end associate

    !$omp end master
    !$omp barrier

  end subroutine Add_SGS_Diffusivity

  !-----------------------------------------------------------------------------
  !> Evaluation on the sigma subgrid-scale model on Cartesian mesh

  subroutine EvalSigmaModel(c_sigma, dx, grad_v, nu)
    real(RNP), intent(in)    :: c_sigma(:)          !< model constant
    real(RNP), intent(in)    :: dx(:)               !< element spacing
    real(RNP), intent(in)    :: grad_v(:,:,:,:,:,:) !< velocity gradient
    real(RNP), intent(inout) :: nu(:,:,:,:,:)       !< diffusivity

    real(RNP) :: delta   ! subgrid characteristic length
    real(RNP) :: nu_sgs  ! subgrid scale viscosity
    integer   :: np, ne
    integer   :: i, j, k, e

    np = size(grad_v,1) ! number of element points per direction
    ne = size(grad_v,4) ! number of elements

    delta   = product(dx / max(np-1,1))**THIRD  ! mean point spacing

    do e = 1, ne
      do k = 1, np
      do j = 1, np
      do i = 1, np
        call SigmaModel(c_sigma(e), delta, grad_v(i,j,k,e,:,:), nu_sgs)
        nu(i,j,k,e,1) = max( nu(i,j,k,e,1) + nu_sgs, ZERO )
        nu(i,j,k,e,2) = max( nu(i,j,k,e,2) + nu_sgs, ZERO )
        nu(i,j,k,e,3) = max( nu(i,j,k,e,3) + nu_sgs, ZERO )
      end do
      end do
      end do
    end do

  end subroutine EvalSigmaModel

  !-----------------------------------------------------------------------------
  !> Evaluation of the model constant in dynamic procedure

  subroutine EvalCSigma(flow_op, u, grad_v, c_sigma_d, d_avg)
    class(FlowOperators), intent(in)    :: flow_op           !< flow operator
    real(RNP),            intent(in)    :: u(:,:,:,:)        !< velocity
    real(RNP),            intent(in)    :: grad_v(:,:,:,:,:) !< velocity gradient
    real(RNP),            intent(out)   :: c_sigma_d         !< dynamic model constant
    real(RNP), optional,  intent(inout) :: d_avg             !< element-averaged dc

    ! following variables for dynamic procedure
    real(RNP), allocatable,save  :: F_c(:,:)            ! filtering matrix
    real(RNP), allocatable       :: uf(:,:,:,:)         ! filtered velocity
    real(RNP), allocatable       :: uc_uc(:,:,:,:)      ! computed u
    real(RNP), allocatable       :: uuf(:,:,:,:)        ! filtered u_u
    real(RNP), allocatable       :: uf_uf(:,:,:,:)
    real(RNP), allocatable       :: dc(:,:,:), df(:,:,:) !computed/filtered D_sigma
    real(RNP), allocatable       :: L(:,:,:,:), M(:,:,:,:)
    real(RNP), allocatable       :: Sc(:,:,:,:), Sf(:,:,:,:)
    real(RNP), allocatable       :: dc_Sf(:,:,:,:), dc_Sc(:,:,:,:)
    real(RNP), allocatable       :: LM(:,:,:), MM(:,:,:)
    real(RNP), allocatable       :: Ms(:) ! weight

    real(RNP)                    :: LM_avg, MM_avg !element averaged terms
    real(RNP)                    :: delta, deltaf  !characteristic lengths
    real(RNP)                    :: c_delta
    integer                      :: po_c ! cut-off degree of the filer
    integer                      :: np
    integer                      :: i,j,k

    ! number of element points per direction, np = po + 1
    np = size(u,1)

    ! cut-off polynomial order of the filter
    po_c = flow_op % control % sgs_dyn_po_c
    if (po_c < 1) then
      po_c = max(1, (np-1)/2)  ! automatic setting to half order
    end if

    associate( mesh  => flow_op % mesh  & ! access to mesh
             , eop   => flow_op % eop_u)  ! access to element operators

    allocate(Ms(np))
    Ms = flow_op%eop_u%w

    allocate(uf, mold = u       )
    allocate(uc_uc(np,np,np,1:6))
    allocate(uf_uf, mold = uc_uc)
    allocate(uuf,   mold = uc_uc)
    allocate(dc(np,np,np)       )
    allocate(df, mold = dc      )
    allocate(L(np,np,np,1:6)    )
    allocate(M,     mold = L    )
    allocate(Sc,    mold = L    )
    allocate(Sf,    mold = L    )
    allocate(dc_Sc, mold = L    )
    allocate(dc_Sf, mold = L    )
    allocate(LM(np,np,np)       )
    allocate(MM, mold = LM      )

    ! step 1: construct the filter matrix ......................................
    if (.not. allocated(F_c)) then
      allocate(F_c(np,np))
      call Filter_Construct(eop, po_c, F_c)
    end if

    ! step 2: component filtering in element level .............................

    ! velocity components
    uf = 0
    do i = 1, 3
      call TPO_AAA_ele(F_c, u(:,:,:,i), uf(:,:,:,i))
    end do

    ! u_i*u_j, u_i(filtered)*u_j(filtered) matrix
    uc_uc = 0
    uf_uf = 0

    ! The matrices uc_uc, uf_uf, uuf, S, S_f, d_S, d_S_f, L, M, LM, MM are
    ! symmetric. Thus only 6 different elements are calculated, with notations:
    ! A(..., 1) = A (...,1,1)
    ! A(..., 2) = A (...,1,2) = A(...,2,1)
    ! A(..., 3) = A (...,1,3) = A(...,3,1)
    ! A(..., 4) = A (...,2,2)
    ! A(..., 5) = A (...,2,3) = A(...,3,2)
    ! A(..., 6) = A (...,3,3)

    uc_uc(:,:,:,1) = u(:,:,:,1)  * u(:,:,:,1)
    uc_uc(:,:,:,2) = u(:,:,:,1)  * u(:,:,:,2)
    uc_uc(:,:,:,3) = u(:,:,:,1)  * u(:,:,:,3)
    uc_uc(:,:,:,4) = u(:,:,:,2)  * u(:,:,:,2)
    uc_uc(:,:,:,5) = u(:,:,:,2)  * u(:,:,:,3)
    uc_uc(:,:,:,6) = u(:,:,:,3)  * u(:,:,:,3)

    uf_uf(:,:,:,1) = uf(:,:,:,1)  * uf(:,:,:,1)
    uf_uf(:,:,:,2) = uf(:,:,:,1)  * uf(:,:,:,2)
    uf_uf(:,:,:,3) = uf(:,:,:,1)  * uf(:,:,:,3)
    uf_uf(:,:,:,4) = uf(:,:,:,2)  * uf(:,:,:,2)
    uf_uf(:,:,:,5) = uf(:,:,:,2)  * uf(:,:,:,3)
    uf_uf(:,:,:,6) = uf(:,:,:,3)  * uf(:,:,:,3)

    ! construct term L and u_i*u_j(filtered)
    uuf = 0
    do i = 1, 6
      call TPO_AAA_ele(F_c, uc_uc(:,:,:,i), uuf(:,:,:,i))
    end do

    ! get dc value for each point
    dc = 0
    do k = 1, np
    do j = 1, np
    do i = 1, np
      call DSigma(grad_v(i,j,k,:,:), dc(i,j,k))
    end do
    end do
    end do

    ! construct Sc, dc*Sc
    Sc    = 0
    dc_Sc = 0
    Sc(:,:,:,1) = grad_v(:,:,:,1,1)
    Sc(:,:,:,2) = HALF* (grad_v(:,:,:,2,1) + grad_v(:,:,:,1,2))
    Sc(:,:,:,3) = HALF* (grad_v(:,:,:,3,1) + grad_v(:,:,:,1,3))
    Sc(:,:,:,4) = grad_v(:,:,:,2,2)
    Sc(:,:,:,5) = HALF* (grad_v(:,:,:,3,2) + grad_v(:,:,:,2,3))
    Sc(:,:,:,6) = grad_v(:,:,:,3,3)

    do i = 1, 6
      dc_Sc(:,:,:,i) = dc(:,:,:) * Sc(:,:,:,i)
    end do

    ! filtered dc
    df = 0
    call TPO_AAA_ele(F_c, dc, df)

    ! Sf, filtered dc*Sc
    Sf    = 0
    dc_Sf = 0

    do i = 1, 6
      call TPO_AAA_ele(F_c, Sc(:,:,:,i), Sf(:,:,:,i))
      call TPO_AAA_ele(F_c, dc_Sc(:,:,:,i), dc_Sf(:,:,:,i))
    end do

    ! step 3: calculate terms L (already done), M ..............................

    !delta and deltaf
    delta  = product(mesh%dx / max(np-1,1))**THIRD
    deltaf = product(mesh%dx / max(po_c,1))**THIRD

    ! L, M
    L = 0
    M = 0
    do i = 1, 6
      L(:,:,:,i) = uuf(:,:,:,i) - uf_uf(:,:,:,i)
      M(:,:,:,i) = (deltaf / delta)**2 * df(:,:,:) * &
                   Sf(:,:,:,i) - dc_Sf(:,:,:,i)
    end do

    ! step 4: sum over space and average over element ..........................
    LM     = 0
    MM     = 0
    LM_avg = 0
    MM_avg = 0
    LM(:,:,:) = L(:,:,:,1) * M(:,:,:,1)     &
              + 2 * L(:,:,:,2) * M(:,:,:,2) &
              + 2 * L(:,:,:,3) * M(:,:,:,3) &
              + L(:,:,:,4) * M(:,:,:,4)     &
              + 2 * L(:,:,:,5) * M(:,:,:,5) &
              + L(:,:,:,6) * M(:,:,:,6)

    MM(:,:,:) = M(:,:,:,1) * M(:,:,:,1)     &
              + 2 * M(:,:,:,2) * M(:,:,:,2) &
              + 2 * M(:,:,:,3) * M(:,:,:,3) &
              + M(:,:,:,4) * M(:,:,:,4)     &
              + 2 * M(:,:,:,5) * M(:,:,:,5) &
              + M(:,:,:,6) * M(:,:,:,6)

    do k = 1, np
    do j = 1, np
    do i = 1, np
      LM_avg = LM_avg + Ms(i)*Ms(j)*Ms(k)*LM(i,j,k)
      MM_avg = MM_avg + Ms(i)*Ms(j)*Ms(k)*MM(i,j,k)
    end do
    end do
    end do
    LM_avg = LM_avg / np**3
    MM_avg = MM_avg / np**3

    c_delta = -HALF*LM_avg / MM_avg
    c_sigma_d = 0
    c_sigma_d = sign(ONE,c_delta) * sqrt(abs(c_delta)) / delta

    ! average of dc for temporal smoothing ....................................

    if (present(d_avg)) then
      d_avg = 0
      do k = 1, np
      do j = 1, np
      do i = 1, np
        d_avg = d_avg + Ms(i)*Ms(j)*Ms(k)*dc(i,j,k)
      end do
      end do
      end do
      d_avg = d_avg / 8
    end if

    ! Finalization .............................................................

    if (allocated(uf)   ) deallocate(uf   )
    if (allocated(uc_uc)) deallocate(uc_uc)
    if (allocated(uuf  )) deallocate(uuf  )
    if (allocated(uf_uf)) deallocate(uf_uf)
    if (allocated(dc)   ) deallocate(dc   )
    if (allocated(df)   ) deallocate(df   )
    if (allocated(L)    ) deallocate(L    )
    if (allocated(M)    ) deallocate(M    )
    if (allocated(Sc)   ) deallocate(Sc   )
    if (allocated(Sf)   ) deallocate(Sf   )
    if (allocated(dc_Sc)) deallocate(dc_Sc)
    if (allocated(dc_Sf)) deallocate(dc_Sf)
    if (allocated(LM)   ) deallocate(LM   )
    if (allocated(MM)   ) deallocate(MM   )
    if (allocated(Ms)   ) deallocate(Ms   )

    end associate

  end subroutine EvalCSigma

  !-----------------------------------------------------------------------------
  !> Get the coefficient from dynamic procedure

  subroutine GetDynamicCoefficient(flow_op, u, grad_v, c_sigma)
    class(FlowOperators), intent(in)  :: flow_op
    real(RNP),            intent(in)  :: u(:,:,:,:,:)
    real(RNP),            intent(in)  :: grad_v(:,:,:,:,:,:)
    real(RNP),            intent(out) :: c_sigma(:)

    ! temporal smoothing variables
    real(RNP), allocatable       :: d_avg(:)
    real(RNP), allocatable, save :: c_sigma_o(:)
    real(RNP), allocatable       :: c2_d_m(:)
    real(RNP), allocatable, save :: c2_d_m_o(:)
    real(RNP), allocatable       :: c_sigma_d(:)
    real(RNP)                    :: eps
    logical                      :: start = .true.

    integer :: i, ne

    ne = size(u,4)

    if (flow_op % control % sgs_smooth) then ! with temporal smoothing .........

      if (.not. allocated(d_avg)    ) allocate (d_avg(ne)    )
      if (.not. allocated(c_sigma_o)) allocate (c_sigma_o(ne))
      if (.not. allocated(c2_d_m)   ) allocate (c2_d_m(ne)   )
      if (.not. allocated(c2_d_m_o) ) allocate (c2_d_m_o(ne) )
      if (.not. allocated(c_sigma_d)) allocate (c_sigma_d(ne))

      eps = flow_op  % control% sgs_smooth_wt

      c_sigma_d = 0
      d_avg = 0
      do i = 1, ne
        call EvalCSigma(  flow_op             &
                        , u(:,:,:,i,:)        &
                        , grad_v(:,:,:,i,:,:) &
                        , c_sigma_d(i)        &
                        , d_avg(i)            )
      end do
      ! start temporal smoothing
      if (start) then
        start = .false.
        c2_d_m = c_sigma_d**2 * d_avg
        c_sigma = c_sigma_d
      else
        c2_d_m  = (1-eps) * c2_d_m_o &
                  + eps * sign(ONE,c_sigma_d) * c_sigma_d**2 * d_avg
        c_sigma = sign(ONE,c2_d_m) * sqrt( abs(c2_d_m / d_avg) )
      end if
      c2_d_m_o = c2_d_m
      c_sigma_o = c_sigma

      deallocate(d_avg)
      deallocate(c2_d_m)
      deallocate(c_sigma_d)

    else ! without temporal smoothing ..........................................

      do i = 1, ne
        call EvalCSigma(  flow_op             &
                        , u(:,:,:,i,:)        &
                        , grad_v(:,:,:,i,:,:) &
                        , c_sigma(i)          )
      end do

    end if

  end subroutine GetDynamicCoefficient

  !-----------------------------------------------------------------------------
  !> build the filter matrix F_c for dynamic procedure

  subroutine Filter_Construct(eop, po_c, F_c)
    use Standard_Element_Operators__1D

    class(StandardElementOperators_1D), intent(in) :: eop
    integer,   intent(in)  :: po_c
    real(RNP), intent(out) :: F_c(:,:)

    real(RNP),allocatable     :: VL(:,:), VL_I(:,:)
    real(RNP),allocatable     :: F(:,:), ws(:,:) ! workspace
    integer                   :: np

    integer :: i,j,k

    np = size(F_c,1)

    allocate(VL,   mold = F_c)
    allocate(VL_I, mold = F_c)
    allocate(F,    mold = F_c)
    allocate(ws,   mold = F_c)

    call eop%Get_Legendre_VDM(VL)
    call eop%Get_Inverse_Legendre_VDM(VL_I)

    ! step 1: construct F ......................................................

    F = 0
    do i = 1,po_c+1
      F(i,i) = 1
    end do

    ! step 2: F_c = VL F VL_I ..................................................
    ws = 0

    do j = 1,np
    do i = 1,np
    do k = 1,np
      ws(i,j) = ws(i,j) + F(i,k)*VL_I(k,j)
    end do
    end do
    end do

    F_c = 0

    do j = 1,np
    do i = 1,np
    do k = 1,np
      F_c(i,j) = F_c(i,j) + VL(i,k)*ws(k,j)
    end do
    end do
    end do

  end subroutine Filter_Construct

  !-----------------------------------------------------------------------------
  !> Tensor product operator 3d AxAxA for dynamic procedure (element-based)

  subroutine TPO_AAA_ele(A, u, v)
    real(RNP), intent(in)  :: A(:,:)   !< 1D operator
    real(RNP), intent(in)  :: u(:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:) !< result

    ! local variables
    real(RNP) :: At( size(A,2), size(A,1) )
    real(RNP) :: z2( size(A,2), size(A,1), size(A,1))
    real(RNP) :: z3( size(A,2), size(A,2), size(A,1))
    real(RNP) :: tmp
    integer   :: na1, na2
    integer   :: i, j, k, p

    ! initialization
    na1 = size(A,1)
    na2 = size(A,2)

    At = transpose(A)

    ! evaluation
    do k = 1, na1
    do j = 1, na2
    do i = 1, na2
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,k) * u(i,j,p)
      end do
      z3(i,j,k) = tmp
    end do
    end do
    end do

    do k = 1, na1
    do j = 1, na1
    do i = 1, na2
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,j) * z3(i,p,k)
      end do
      z2(i,j,k) = tmp
    end do
    end do
    end do

    do k = 1, na1
    do j = 1, na1
    do i = 1, na1
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,i) * z2(p,j,k)
      end do
      v(i,j,k) = tmp
    end do
    end do
    end do

  end subroutine TPO_AAA_ele

  !=============================================================================

end module CART__ISP_Flow__SGS_Diffusivity
