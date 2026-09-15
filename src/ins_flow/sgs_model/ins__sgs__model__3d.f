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

!> summary:  Subgrid scale model for incompressible flow
!> author:   Yang Liu, Moritz Kreuseler, Joerg Stiller
!> date:     2026/09/06
!===============================================================================

module INS__SGS__Model__3D
  use Kind_Parameters, only: RNP
  use Constants
  use XMPI
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Spectral_Element_Mesh__3D
  use Smooth_Mesh_Data__3D

  implicit none
  private

  public :: INS_SGS_Options_3D
  public :: INS_SGS_Model_3D

  !-----------------------------------------------------------------------------
  !> SGS options

  type INS_SGS_Options_3D
    integer   :: model     = 0      !< SGS model, 0/1/2: none/Smagorinsky/Sigma
    integer   :: length    = 1      !< element length scale, 1/2/3: mean/max/min
    integer   :: smooth    = 2      !< smoothing, -1/0/1/2: none/avg/linear/half
    logical   :: dynamic   = .true. !< switch to dynamic model
    real(RNP) :: c_static  = 1.35   !< coefficient of static model
  contains
    procedure :: Bcast => Bcast_INS_SGS_Options_3D
  end type INS_SGS_Options_3D

  !-----------------------------------------------------------------------------
  !> SGS model

  type INS_SGS_Model_3D
    integer   :: model    !< SGS model, 0/1/2: none/Smagorinsky/Sigma
    integer   :: length   !< element length scale, 1/2/3: mean/max/min
    integer   :: smooth   !< smoothing, -1/0/1/2: none/avg/linear/half degree
    logical   :: dynamic  !< switch to dynamic model
    real(RNP) :: c_static !< coefficient of static model
  contains
    procedure :: Get_SGS_Viscosity
  end type INS_SGS_Model_3D

  ! constructor interface
  interface INS_SGS_Model_3D
    module procedure New_INS_SGS_Model_3D
  end interface

  real(RNP), parameter :: EPS = epsilon(ONE)

contains

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of SGS optiones

  subroutine Bcast_INS_SGS_Options_3D(this, root, comm)
    class(INS_SGS_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % model   , root, comm)
    call XMPI_Bcast(this % length  , root, comm)
    call XMPI_Bcast(this % smooth  , root, comm)
    call XMPI_Bcast(this % dynamic , root, comm)
    call XMPI_Bcast(this % c_static, root, comm)

  end subroutine Bcast_INS_SGS_Options_3D

  !-----------------------------------------------------------------------------
  !> Returns a new SGS model

  type(INS_SGS_Model_3D) function New_INS_SGS_Model_3D(opt) result(this)
    class(INS_SGS_Options_3D), intent(in) :: opt

    call Init_INS_SGS_Model_3D(this, opt)

  end function New_INS_SGS_Model_3D

  !-----------------------------------------------------------------------------
  !> Initialization of SGS model

  subroutine Init_INS_SGS_Model_3D(this, opt)
    class(INS_SGS_Model_3D), intent(inout) :: this
    class(INS_SGS_Options_3D), intent(in) :: opt

    this % model    = opt % model
    this % length   = opt % length
    this % smooth   = opt % smooth
    this % dynamic  = opt % dynamic
    this % c_static = opt % c_static

  end subroutine Init_INS_SGS_Model_3D

  !-----------------------------------------------------------------------------
  !> Computation of SGS viscosity

  subroutine Get_SGS_Viscosity(this, sem, v, nu)
    class(INS_SGS_Model_3D), intent(in) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem  !< spectral element mesh
    real(RNP), contiguous, intent(in)  :: v(:,:,:,:,:)!< flow velocity
    real(RNP), contiguous, intent(out) :: nu(:,:,:,:) !< SGS viscosity

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    real(RNP), allocatable, save :: vp(:,:,:,:,:)       ! velocity traces v⁺
    real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:) ! velocity gradient ∇v
    real(RNP), allocatable, save :: A(:,:)              ! test filter matrix

    real(RNP), allocatable :: D(:,:,:) ! model operator

    real(RNP) :: dx(3), delta, c_delta2, rr
    integer   :: ne, np, pf
    integer   :: e, i, j, k

    associate( mesh => sem % mesh   &
             , eop  => sem % std_op )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      np = size(v, 1)
      ne = mesh % n_elem

      if (this % dynamic) then
        pf = ceiling(HALF * eop%po)
        rr = (eop%po / max(pf, 1))**2
      end if

      allocate(D(np,np,np))

      !$omp master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      allocate(grad_v(np,np,np,ne,3,3))
      allocate(vp(np,np,6,ne,3))
      allocate(A(np,np))
      if (this % dynamic) then
        call eop % Get_Legendre_CutoffFilter(pf, A)
      end if
      !$omp end master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !$omp barrier

      ! velocity gradient ::::::::::::::::::::::::::::::::::::::::::::::::::::::

      call GetOuterVectorTraces_3D(mesh, v(:,:,:,:,1:3), vp)

      do k = 1, 3
        call TPO_Grad(eop, sem, v(:,:,:,:,k), vp(:,:,:,:,k), grad_v(:,:,:,:,:,k))
      end do

      ! SGS viscosity ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp do
      do e = 1, ne

        ! filter width .........................................................

        call mesh % element(e) % GetCuboidDimensions(dx)

        select case(this % length)
        case(2)
          delta = maxval(dx)
        case(3)
          delta = minval(dx)
        case default
          delta = ( dx(1)**2 + dx(2)**2 + dx(3)**2 ) ** THIRD
        end select
        delta = delta / max(eop % po, 1)

        ! model operator .......................................................

        select case(this % model)
        case(1)
          do concurrent(i = 1:np, j = 1:np, k = 1:np)
            D(i,j,k) = D_Smagorinsky(grad_v(i,j,k,e,:,:))
          end do
        case(2)
          do concurrent(i = 1:np, j = 1:np, k = 1:np)
            D(i,j,k) = D_Sigma(grad_v(i,j,k,e,:,:))
          end do
        end select

        ! model coefficient (CΔ)² ..............................................

        if (this % dynamic) then
          call DynamicProcedure(e, eop%w, rr, A, D, v, grad_v, c_delta2)
        else
          c_delta2 = (this % c_static * delta)**2
        end if

        ! SGS viscosity ........................................................

        nu(:,:,:,e) = c_delta2 * max(D, ZERO)

      end do

      ! smoothing ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      block
        integer :: filter
        integer :: order

        select case(this % smooth)
        case(0)
          filter = 0
          order  = 0
        case(1)
          filter = 1
          order  = 1
        case(2)
          filter = 1
          order  = ceiling(HALF * eop%po)
        end select

        if (this % smooth >= 0) then
          call SmoothMeshData_3D(mesh, eop, nu, filter, order)
        end if

      end block

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(grad_v, vp, A)
      !$omp end master

    end associate

  end subroutine Get_SGS_Viscosity

  !-----------------------------------------------------------------------------
  !> Differential operator of the Smagorinsky model

  pure real(RNP) function D_Smagorinsky(grad_v) result(D)
    real(RNP), intent(in) :: grad_v(:,:) ! velocity gradient (3,3)
    integer :: i, j

    D = ZERO
    do j = 1, 3
    do i = 1, 3
      D = D + (grad_v(i,j) + grad_v(j,i))**2
    end do
    end do
    D = sqrt(HALF * D)

  end function D_Smagorinsky

  !-----------------------------------------------------------------------------
  !> Differential operator of the Sigma model

  pure real(RNP) function D_Sigma(grad_v) result(D)
    real(RNP), intent(in) :: grad_v(:,:) ! velocity gradient (3,3)

    real(RNP) :: G(3,3), GG(3,3)
    real(RNP) :: alpha_1, sigma_1, I_1
    real(RNP) :: alpha_2, sigma_2, I_2
    real(RNP) :: alpha_3, sigma_3, I_3
    real(RNP) :: a, b, c
    integer   :: i, j, k

    ! initialization ...........................................................

    ! G = ∇vᵀ·∇v
    G = 0
    do j = 1, 3
    do i = 1, 3
    do k = 1, 3
      G(i,j) = G(i,j) + grad_v(k,i) * grad_v(k,j)
    end do
    end do
    end do

    ! GG = G²
    GG = 0
    do j = 1, 3
    do i = 1, 3
    do k = 1, 3
      GG(i,j) = GG(i,j) + G(i,k) * G(k,j)
    end do
    end do
    end do

    ! invariants ...............................................................

    ! I₁ = Tr(G)
    I_1 = G(1,1) + G(2,2) + G(3,3)

    ! I₂ = (Tr(G)² - Tr(G²)) / 2
    I_2 =(I_1**2 - (GG(1,1) + GG(2,2) + GG(3,3))) * HALF

    ! I₃ = det(G)
    I_3 = G(1,1) * G(2,2) * G(3,3)  &
        + G(1,2) * G(2,3) * G(3,1)  &
        + G(1,3) * G(2,1) * G(3,2)  &
        - G(3,1) * G(2,2) * G(1,3)  &
        - G(1,1) * G(2,3) * G(3,2)  &
        - G(1,2) * G(2,1) * G(3,3)

    ! angles ...................................................................

    ! α₁ = (I₁/3)² - I₂/3
    a = I_1 * THIRD
    alpha_1 = a**2 - I_2 * THIRD

    ! α₂ = (I₁/3)³ - I₁I₂/6 + I₃/2
    alpha_2 = a**3 + (-a * I_2 + I_3) * HALF

    ! α₃ = arccos(α₂ / (√α₁)³) / 3
    b = sqrt(max(alpha_1, ZERO))
    c = alpha_2 / max(b**3, EPS)
    c = sign(ONE, c) * min(ONE, abs(c))  ! ensure c ∈ [-1,1]
    alpha_3 = acos(c) * THIRD

    ! singular values ..........................................................

    c = PI * THIRD

    sigma_1 = sqrt( max(ZERO, a + 2*b * cos(    alpha_3)) )
    sigma_2 = sqrt( max(ZERO, a - 2*b * cos(c + alpha_3)) )
    sigma_3 = sqrt( max(ZERO, a - 2*b * cos(c - alpha_3)) )

    ! differential operator ....................................................

    D = sigma_3 * (sigma_1 - sigma_2) * (sigma_2 - sigma_3) &
      / max(sigma_1**2, EPS)

  end function D_Sigma

  !-----------------------------------------------------------------------------
  !> Dynamic evaluation of the SGS model coefficient

  pure subroutine DynamicProcedure(e, w, rr, A, D, v, grad_v, c_delta2)
    integer,   intent(in)  :: e                   !< element ID
    real(RNP), intent(in)  :: w(:)                !< 1D quadrature weights
    real(RNP), intent(in)  :: rr                  !< squared filter ratio
    real(RNP), intent(in)  :: A(:,:)              !< test filter matrix
    real(RNP), intent(in)  :: D(:,:,:)            !< model operator
    real(RNP), intent(in)  :: v(:,:,:,:,:)        !< velocity v
    real(RNP), intent(in)  :: grad_v(:,:,:,:,:,:) !< velocity gradient ∇v
    real(RNP), intent(out) :: c_delta2            !< dynamic coefficient

    contiguous :: w, A, D, v, grad_v

    ! internal variables .......................................................

    ! current
    real(RNP), dimension(size(w), size(w), size(w), 6) :: vv, DS, L, M, S

    ! filtered
    real(RNP), dimension(size(w), size(w), size(w))    :: Df
    real(RNP), dimension(size(w), size(w), size(w), 3) :: vf
    real(RNP), dimension(size(w), size(w), size(w), 6) :: vvf, DSf, Sf

    real(RNP) :: LM, MM

    integer :: i, j, k, np

    ! initialization ...........................................................

    np = size(w)

    ! filtered velocity
    do k = 1, 3
      call TPO_AAA(A, v(:,:,:,e,k), vf(:,:,:,k))
    end do

    ! momentum flux tensor, skipping symmetric components
    vv(:,:,:,1) = v(:,:,:,e,1) * v(:,:,:,e,1)
    vv(:,:,:,2) = v(:,:,:,e,1) * v(:,:,:,e,2)
    vv(:,:,:,3) = v(:,:,:,e,1) * v(:,:,:,e,3)
    vv(:,:,:,4) = v(:,:,:,e,2) * v(:,:,:,e,2)
    vv(:,:,:,5) = v(:,:,:,e,2) * v(:,:,:,e,3)
    vv(:,:,:,6) = v(:,:,:,e,3) * v(:,:,:,e,3)

    ! filtered momentum flux tensor
    do k = 1, 6
      call TPO_AAA(A, vv(:,:,:,k), vvf(:,:,:,k))
    end do

    ! deformation tensor, skipping symmetric components
    S(:,:,:,1) = grad_v(:,:,:,e,1,1)                                 ! S₁₁
    S(:,:,:,2) = HALF * (grad_v(:,:,:,e,1,2) + grad_v(:,:,:,e,2,1))  ! S₁₂
    S(:,:,:,3) = HALF * (grad_v(:,:,:,e,1,3) + grad_v(:,:,:,e,3,1))  ! S₁₃
    S(:,:,:,4) = grad_v(:,:,:,e,2,2)                                 ! S₂₂
    S(:,:,:,5) = HALF * (grad_v(:,:,:,e,2,3) + grad_v(:,:,:,e,3,2))  ! S₂₃
    S(:,:,:,6) = grad_v(:,:,:,e,3,3)                                 ! S₃₃

    do k = 1, 6
      ! filtered deformation tensor
      call TPO_AAA(A, S(:,:,:,k), Sf(:,:,:,k))
      ! scaled deformation tensor
      DS(:,:,:,k) = D * S(:,:,:,k)
      ! filtered scaled deformation tensor
      call TPO_AAA(A, DS(:,:,:,k), DSf(:,:,:,k))
    end do

    ! filtered model operator
    call TPO_AAA(A, D, Df)

    ! test operators ..........................................................

    ! modified Leonard stress tensor, skipping symmetric components
    L(:,:,:,1) = vvf(:,:,:,1) - vf(:,:,:,1) * vf(:,:,:,1)  ! L₁₁
    L(:,:,:,2) = vvf(:,:,:,2) - vf(:,:,:,1) * vf(:,:,:,2)  ! L₁₂
    L(:,:,:,3) = vvf(:,:,:,3) - vf(:,:,:,1) * vf(:,:,:,3)  ! L₁₃
    L(:,:,:,4) = vvf(:,:,:,4) - vf(:,:,:,2) * vf(:,:,:,2)  ! L₂₂
    L(:,:,:,5) = vvf(:,:,:,5) - vf(:,:,:,2) * vf(:,:,:,3)  ! L₂₃
    L(:,:,:,6) = vvf(:,:,:,6) - vf(:,:,:,3) * vf(:,:,:,3)  ! L₃₃

    ! scaled SGS stress in the test window
    do k = 1, 6
      M(:,:,:,k) = DSf(:,:,:,k) - rr * Df(:,:,:) * Sf(:,:,:,k)
    end do

    ! model coefficient ........................................................

    LM = 0
    MM = 0
    do k = 1, np
    do j = 1, np
    do i = 1, np

      LM = LM + w(i)*w(j)*w(k) * ( L(i,j,k,1) * M(i,j,k,1)      &
                                 + L(i,j,k,2) * M(i,j,k,2) * 2  &
                                 + L(i,j,k,3) * M(i,j,k,3) * 2  &
                                 + L(i,j,k,4) * M(i,j,k,4)      &
                                 + L(i,j,k,5) * M(i,j,k,5) * 2  &
                                 + L(i,j,k,6) * M(i,j,k,6)      )

      MM = MM + w(i)*w(j)*w(k) * ( M(i,j,k,1) * M(i,j,k,1)      &
                                 + M(i,j,k,2) * M(i,j,k,2) * 2  &
                                 + M(i,j,k,3) * M(i,j,k,3) * 2  &
                                 + M(i,j,k,4) * M(i,j,k,4)      &
                                 + M(i,j,k,5) * M(i,j,k,5) * 2  &
                                 + M(i,j,k,6) * M(i,j,k,6)      )
    end do
    end do
    end do

    c_delta2 = HALF * max(LM, ZERO) / max(MM, EPS)

  end subroutine DynamicProcedure

  !----------------------------------------------------------------------------
  !> Tensor product operator AxAxA for dynamic procedure

  pure subroutine TPO_AAA(A, u, v)
    real(RNP), intent(in)  :: A(:,:)   !< 1D operator
    real(RNP), intent(in)  :: u(:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:) !< result

    ! internal variables ......................................................

    real(RNP) :: At( size(A,2), size(A,1) )
    real(RNP) :: z2( size(A,2), size(A,1), size(A,1))
    real(RNP) :: z3( size(A,2), size(A,2), size(A,1))
    real(RNP) :: tmp
    integer   :: na1, na2
    integer   :: i, j, k, p

    ! initialization ..........................................................

    na1 = size(A,1)
    na2 = size(A,2)

    At = transpose(A)

    ! evaluation ..............................................................

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

  end subroutine TPO_AAA

  !=============================================================================

end module INS__SGS__Model__3D
