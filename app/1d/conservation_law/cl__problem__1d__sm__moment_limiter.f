submodule(CL__Problem__1D) SM_Moment_Limiter
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of moment limiter

  module subroutine MomentLimiter(this, cl_operator, u)
    class(CL_Problem_1D),  intent(in)    :: this
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP), contiguous, intent(inout) :: u(0:,:,:)

    real(RNP), allocatable, save :: u_o(:,:,:)

    real(RNP), allocatable :: VL(:,:), VL_inv(:,:)
    real(RNP) :: s
    integer   :: e, i

    associate( nc => this % nc                    &
             , ne => cl_operator % ne             &
             , po => cl_operator % eop % po       &
             , activity => cl_operator % activity )

      ! preliminaries ..........................................................

      allocate(u_o(0:po,0:ne+1,nc))

      allocate(VL(0:po,0:po), VL_inv(0:po,0:po))
      call cl_operator % eop % Get_Legendre_VDM(VL)
      call cl_operator % eop % Get_Inverse_Legendre_VDM(VL_inv)

      ! transform to Legendre basis ............................................

      do e = 1, ne
        if (activity(e) < 0) cycle
        u  (0:po,e,1:nc) = matmul(VL_inv, u(0:po,e,1:nc))
        u_o(0:po,e,1:nc) = u(0:po,e,1:nc)
      end do

      ! left boundary
      select case(this % bc(1))
      case('P')
        ! use periodicity
        u_o(0:po,0,1:nc) = u_o(0:po,ne,1:nc)
      case default
        ! apply reflection
        s = 1
        do i = 0, po
          u_o(i,0,1:nc) = s * u_o(i,1,1:nc)
          s = -s
        end do
      end select

      ! right boundary
      select case(this % bc(2))
      case('P')
        ! use periodicity
        u_o(0:po,ne+1,1:nc) = u_o(0:po,1,1:nc)
      case default
        s = 1
        do i = 0, po
          u_o(i,ne+1,1:nc) = s * u_o(i,ne,1:nc)
          s = -s
        end do
      end select

      ! apply limiter ..........................................................

      do e = 1, ne
        if (activity(e) < 1) cycle
        if (nc == 1) then
          call BurbeauLimiter_Scalar(u(:,e,1), u_o(:,e-1,1), u_o(:,e+1,1))
        else
          call BurbeauLimiter_System(this, u(:,e,:), u_o(:,e-1,:), u_o(:,e+1,:))
        end if
      end do

      ! transform to nodal basis ...............................................

      do e = 1, ne
        if (activity(e) < 1) cycle
        u(0:po,e,1:nc) = matmul(VL, u(0:po,e,1:nc))
      end do

      ! finalization ...........................................................

      deallocate(u_o)

    end associate

  end subroutine MomentLimiter

  !-----------------------------------------------------------------------------
  !> Slope limiter of Burbeau et al. 2001 for scalar conservation laws
  !>
  !> Source: <https://doi.org/10.1006/jcph.2001.6718>

  subroutine BurbeauLimiter_Scalar(u, u_l, u_r)
    real(RNP), intent(inout) :: u  (0:) !< element solution being limited
    real(RNP), intent(in)    :: u_l(0:) !< unlimited solution in left element
    real(RNP), intent(in)    :: u_r(0:) !< unlimited solution in right element

    real(RNP) :: du_l, du_r, du_0, du_m, du_max
    real(RNP) :: c
    integer   :: m, po

    po = ubound(u,1)

    do m = po-1, 0, -1
      c = ONE / (2*m + 1)
      du_l = c * (u  (m) - u_l(m))
      du_r = c * (u_r(m) - u  (m))
      du_0 = u(m+1)
      du_m = minmod(du_0, du_r, du_l)
      if (du_m == du_0) exit
      du_r = du_r - u_r(m+1)
      du_l = du_l - u_l(m+1)
      du_max = minmod(du_0, du_r, du_l)
      u(m+1) = maxmod(du_m, du_max)
    end do

  end subroutine BurbeauLimiter_Scalar

  !-----------------------------------------------------------------------------
  !> Slope limiter of Burbeau et al. 2001 for conservation systems
  !>
  !> Source: <https://doi.org/10.1006/jcph.2001.6718>

  subroutine BurbeauLimiter_System(this, u, u_l, u_r)
    class(CL_Problem_1D), intent(in) :: this
    real(RNP), intent(inout) :: u  (0:,:) !< element solution being limited
    real(RNP), intent(in)    :: u_l(0:,:) !< unlimited solution in left element
    real(RNP), intent(in)    :: u_r(0:,:) !< unlimited solution in right element

    real(RNP), dimension(this%nc,this%nc) :: L, R
    real(RNP), dimension(this%nc) :: du_l, du_r, dz_l, dz_r, dz_0, dz_m, dz_max
    logical  , dimension(this%nc) :: skip

    real(RNP) :: c
    integer   :: m, nc, po

    skip = .false.

    po = ubound(u,1)
    nc = this % nc

    call this % ConvectiveEigensystem(u(0,1:nc), R=R, L=L)

    do m = po-1, 0, -1
      c = ONE / (2*m + 1)
      du_l = c * (u  (m,1:nc) - u_l(m,1:nc))
      du_r = c * (u_r(m,1:nc) - u  (m,1:nc))
      dz_r = matmul(L, du_r)
      dz_l = matmul(L, du_l)
      dz_0 = matmul(L, u(m+1,1:nc))
      dz_m = minmod(dz_0, dz_r, dz_l)
      skip = skip .or. dz_m == dz_0
      if (all(skip)) exit
      dz_r = dz_r - matmul(L, u_r(m+1,1:nc))
      dz_l = dz_l - matmul(L, u_l(m+1,1:nc))
      where (skip)
        dz_max = dz_0
      elsewhere
        dz_max = minmod(dz_0, dz_r, dz_l)
      end where
      u(m+1,1:nc) = matmul(R, maxmod(dz_m, dz_max))
    end do

  end subroutine BurbeauLimiter_System

  !-----------------------------------------------------------------------------
  ! MinMod function

  elemental function MinMod(a_1, a_2, a_3) result(a)
    real(RNP), intent(in) :: a_1
    real(RNP), intent(in) :: a_2
    real(RNP), intent(in) :: a_3
    real(RNP) :: a

    real(RNP) :: a_min, a_max

    a_min = min(a_1, a_2, a_3)
    a_max = max(a_1, a_2, a_3)

    if (a_max <= 0) then
      a = a_max
    else if (a_min >= 0) then
      a = a_min
    else
      a = 0
    end if

  end function MinMod

  !-----------------------------------------------------------------------------
  ! MaxMod function

  elemental function MaxMod(a_1, a_2) result(a)
    real(RNP), intent(in) :: a_1
    real(RNP), intent(in) :: a_2
    real(RNP) :: a

    real(RNP) :: a_min, a_max

    a_min = min(a_1, a_2)
    a_max = max(a_1, a_2)

    if (a_max <= 0) then
      a = a_min
    else if (a_min >= 0) then
      a = a_max
    else
      a = 0
    end if

  end function MaxMod

  !=============================================================================

end submodule SM_Moment_Limiter
