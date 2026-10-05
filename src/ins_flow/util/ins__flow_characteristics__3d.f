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

!> summary:  Evaluation of flow characteristics
!> author:   Joerg Stiller
!> date:     2023/02/17
!===============================================================================

module INS__Flow_Characteristics__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use Array_Reductions
  use XMPI

  use Trace_Operators__3D
  use TPO__Div__3D
  use Boundary_Variable__3D
  use Volume_Integrals__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_FlowCharacteristics_3D

  !-----------------------------------------------------------------------------
  !> Incompressible flow statistics

  type INS_FlowCharacteristics_3D

    real(RNP) :: t      !< problem time
    real(RNP) :: dt     !< time step width
    real(RNP) :: dx_min !< minimum cuboid extension
    real(RNP) :: dx_max !< maximum cuboid extension
    integer   :: po     !< polynomial order
    real(RNP) :: err_v  !< L² velocity error
    real(RNP) :: err_p  !< L² velocity error
    real(RNP) :: div_v  !< L² velocity divergence
    real(RNP) :: v_max  !< maximum velocity
    real(RNP) :: e_kin  !< kinetic energy per unit mass

    ! only if dissipation requested (diss = T):
    real(RNP) :: phi_c  !< convection
    real(RNP) :: phi_d  !< diffusion with μ and ν as given
    real(RNP) :: phi_ds !< diffusion with μ=0 and ν as given
    real(RNP) :: phi_d0 !< diffusion with μ=0 and ν=ν₀
    real(RNP) :: phi_s  !< body force

    integer, private :: part = -1
    logical, private :: normalized      = .false. !< normalized by volume
    logical, private :: has_errors      = .false. !< errors computed
    logical, private :: has_dissipation = .false. !< dissipation computed

  contains
    procedure :: Evaluate
    procedure :: PrintHeader
    procedure :: PrintValues
  end type INS_FlowCharacteristics_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of incompressible flow characteristics

  subroutine Evaluate(this, ins_op, t, mu, nu, u, dt, volume, diss, leaf)
    class(INS_FlowCharacteristics_3D), intent(inout) :: this
    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP),              intent(in) :: t
    real(RNP), contiguous,  intent(in) :: mu(:,:,:,:)
    real(RNP), contiguous,  intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous,  intent(in) :: u(:,:,:,:,:)
    real(RNP), optional,    intent(in) :: dt
    real(RNP), optional,    intent(in) :: volume
    logical,   optional,    intent(in) :: diss   !< compute disspation [F]
    logical,   optional,    intent(in) :: leaf   !< only leaf elements [F]

    type(BoundaryVariable_3D), allocatable, save :: bv_x(:), bv_u(:)
    real(RNP), allocatable, save :: vp(:,:,:,:,:)
    real(RNP), allocatable, save :: sp(:,:,:,:,:)
    real(RNP), allocatable, save :: mm(:,:,:,:)
    real(RNP), allocatable, save :: s0(:,:,:,:)
    real(RNP), allocatable, save :: s1(:,:,:,:)
    real(RNP), allocatable, save :: w(:,:,:,:,:)
    real(RNP), save :: dx_min_loc = huge(ONE)
    real(RNP), save :: dx_max_loc = 0
    real(RNP), save :: v_max      = 0
    real(RNP), save :: vv_max     = 0
    real(RNP), save :: e_kin      = 0
    real(RNP), save :: phi_c      = 0
    real(RNP), save :: phi_d      = 0
    real(RNP), save :: phi_ds     = 0
    real(RNP), save :: phi_d0     = 0
    real(RNP), save :: phi_s      = 0

    real(RNP) :: div_v, err_v, err_p, vv
    real(RNP) :: dx(3)

    logical :: complete
    integer :: b, e, i, j, k, na, nc, ne, np

    if (present(leaf)) then
      complete = .not. leaf
    else
      complete = .true.
    end if

    associate( problem => ins_op % problem             &
             , mesh    => ins_op % mesh                &
             , eop     => ins_op % eop_u               &
             , x       => ins_op % sem_u % metrics % x )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (mesh % part < 0) then
        !$omp barrier
        return
      end if

      np = size(u,1)
      ne = size(u,4)
      nc = size(u,5)
      na = mesh % n_elem_active

      !$omp master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      this % t          = t
      this % part       = mesh % part
      this % normalized = present(volume)
      this % has_errors = problem % HasExactSolution()

      if (present(diss)) then
        this % has_dissipation = diss
      else
        this % has_dissipation = .false.
      end if

      if (present(dt)) then
        this % dt = dt
      else
        this % dt = -1
      end if

      this % po    = eop % po
      this % err_v = -1
      this % err_p = -1

      allocate(vp(np,np,6,ne,3))
      allocate(sp(np,np,6,ne,3))
      allocate(mm(np,np,np,ne)  , source = ZERO )
      allocate(s0(np,np,np,ne)  , source = ZERO )
      allocate(s1(np,np,np,ne)  , source = ONE )
      allocate(w(np,np,np,ne,nc), source = ZERO)

      if (this % has_dissipation) then
        allocate(bv_x (mesh % n_bound) )
        allocate(bv_u (mesh % n_bound) )
        do b = 1, mesh % n_bound
          call bv_u(b) % Init(mesh % boundary(b), eop % po, nc = 4)
          call bv_x(b) % Init(mesh % boundary(b), eop % po, nc = 3)
          call bv_x(b) % Extract(x)
        end do
      end if

      !$omp end master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      ! mesh spacing :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp do private(dx) reduction(min:dx_min_loc) reduction(max:dx_max_loc)
      do e = 1, na
        if (complete .or. mesh % element(e) % IsLeaf()) then
          call mesh % element(e) % GetCuboidDimensions(dx)
          dx_min_loc = min(dx_min_loc, dx(1), dx(2), dx(3))
          dx_max_loc = max(dx_max_loc, dx(1), dx(2), dx(3))
        end if
      end do

      !$omp master
      call XMPI_Allreduce(dx_min_loc, this%dx_min, MPI_MIN, mesh%comm_parts)
      call XMPI_Allreduce(dx_max_loc, this%dx_max, MPI_MAX, mesh%comm_parts)
      !$omp end master

      ! errors :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (this % has_errors) then
        associate(q => w(:,:,:,:,4))

          call problem % GetExactSolution(x, t, w)         !
          call MergeArrays(-ONE, w, ONE, u, multi=.true.)  ! w  = u  - u_ex
          if (complete .and. ne == na) then
            ! remove mean pressure error
            call CalibrateArray(q, comm = mesh%comm_parts)
          end if

          ! pressure ...........................................................

          !$omp do
          do e = 1, na
            q(:,:,:,e) = q(:,:,:,e)**2
          end do
          call GetVolumeIntegral(ins_op%sem_u, q, err_p, leaf)
          err_p = sqrt(err_p)

          ! velocity ...........................................................

          !$omp do
          do e = 1, na
            q(:,:,:,e) = w(:,:,:,e,1)**2 + w(:,:,:,e,2)**2 + w(:,:,:,e,3)**2
          end do
          call GetVolumeIntegral(ins_op%sem_u, q, err_v, leaf)
          err_v = sqrt(err_v)

        end associate
      end if

      ! divergence :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      associate(v => u(:,:,:,:,1:3), q => w(:,:,:,:,4))
        call GetOuterVectorTraces_3D(mesh, v, vp)   ! vp = v⁺
        call TPO_Div(eop, ins_op%sem_u, v, vp, q)   ! q = ∇⋅v
        !$omp do
        do e = 1, na
          q(:,:,:,e) = q(:,:,:,e) ** 2
        end do
        !$omp end do nowait
        !$omp do
        do e = na+1, ne
          q(:,:,:,e) = ZERO
        end do
        call GetVolumeIntegral(ins_op%sem_u, q, div_v, leaf)
        div_v = sqrt(div_v)
      end associate

      ! maximum velocity and energy ::::::::::::::::::::::::::::::::::::::::::::

      associate(v => u(:,:,:,:,1:3), q => w(:,:,:,:,4))

        !$omp do reduction(max:vv_max)
        do e = 1, na
          if (complete .or. mesh % element(e) % IsLeaf()) then
            do k = 1, np
            do j = 1, np
            do i = 1, np
              vv = u(i,j,k,e,1)**2 + u(i,j,k,e,2)**2 + u(i,j,k,e,3)**2
              vv_max = max(vv_max, vv)
              q(i,j,k,e) = vv
            end do
            end do
            end do
          else
            q(:,:,:,e) = ZERO
          end if
        end do

        !$omp master
        call XMPI_Allreduce(sqrt(vv_max), v_max, MPI_MAX, mesh%comm_parts)
        !$omp end master

        call GetVolumeIntegral(ins_op%sem_u, q, e_kin, leaf)

      end associate

      ! dissipation ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (this % has_dissipation) then
        associate(v => u(:,:,:,:,1:3), q => w(:,:,:,:,4))

          ! BC .................................................................

          do b = 1, mesh % n_bound
            select case(problem % bc_v(b))
            case('D')
              call problem % GetBoundaryValues(b, bv_x(b)%val, t, bv_u(b)%val)
            end select
          end do

          ! diffusion with μ and ν as given ....................................

          call ins_op % GetDiffusionTerm( mu, nu, bv_u, w, v, vp, sp &
                                        , xout = .true. )

          !$omp do reduction(+:phi_d)
          do e = 1, na
            if (complete .or. mesh % element(e) % IsLeaf()) then
              phi_d = phi_d + sum(v(:,:,:,e,1) * w(:,:,:,e,1)) &
                            + sum(v(:,:,:,e,2) * w(:,:,:,e,2)) &
                            + sum(v(:,:,:,e,3) * w(:,:,:,e,3))
            end if
          end do

          ! diffusion with μ = 0 and ν as given (shear only) ...................

          call ins_op % GetDiffusionTerm( s0, nu, bv_u, w, v, vp, sp &
                                        , xout = .true. )

          !$omp do reduction(+:phi_ds)
          do e = 1, na
            if (complete .or. mesh % element(e) % IsLeaf()) then
              phi_ds = phi_ds + sum(v(:,:,:,e,1) * w(:,:,:,e,1)) &
                              + sum(v(:,:,:,e,2) * w(:,:,:,e,2)) &
                              + sum(v(:,:,:,e,3) * w(:,:,:,e,3))
            end if
          end do

          ! diffusion with μ = 0 and ν = ν₀ (constant part) ....................

          if (ins_op % nu_0 > 0) then

            call ins_op % GetDiffusionTerm( s0, s1, bv_u, w, v, vp, sp &
                                          , xout = .true. )

            !$omp do reduction(+:phi_d0)
            do e = 1, na
              if (complete .or. mesh % element(e) % IsLeaf()) then
                phi_d0 = phi_d0 + sum(v(:,:,:,e,1) * w(:,:,:,e,1)) &
                                + sum(v(:,:,:,e,2) * w(:,:,:,e,2)) &
                                + sum(v(:,:,:,e,3) * w(:,:,:,e,3))
              end if
            end do
            !$omp master
            phi_d0 = phi_d0 * ins_op % nu_0
            !$omp end master

          end if

          ! convection .........................................................

          call ins_op % GetConvectionTerm(v, vp, w)

          !$omp do reduction(+:phi_c)
          do e = 1, na
            if (complete .or. mesh % element(e) % IsLeaf()) then
              phi_c = phi_c + sum(v(:,:,:,e,1) * w(:,:,:,e,1)) &
                            + sum(v(:,:,:,e,2) * w(:,:,:,e,2)) &
                            + sum(v(:,:,:,e,3) * w(:,:,:,e,3))
            end if
          end do

          ! body force .........................................................

          call problem % GetExternalSources(x, t, w)
          call ins_op % sem_u % Get_DG_DiagonalMassMatrix(mm)

          !$omp do reduction(+:phi_c)
          do e = 1, na
            if (complete .or. mesh % element(e) % IsLeaf()) then
              phi_s = phi_s + sum(mm(:,:,:,e) * v(:,:,:,e,1) * w(:,:,:,e,1)) &
                            + sum(mm(:,:,:,e) * v(:,:,:,e,2) * w(:,:,:,e,2)) &
                            + sum(mm(:,:,:,e) * v(:,:,:,e,3) * w(:,:,:,e,3))
            end if
          end do
        end associate
      end if

      !$omp master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      ! completion .............................................................

      if (this % normalized) then
        this % err_v = err_v / sqrt(volume)
        this % err_p = err_p / sqrt(volume)
        this % div_v = div_v / sqrt(volume)
        this % e_kin = e_kin / (2 * volume)
      else
        this % err_v = err_v
        this % err_p = err_p
        this % div_v = div_v
        this % e_kin = e_kin / 2
      end if
      this % v_max = v_max

      if (this % has_dissipation) then
        call XMPI_Allreduce(phi_c , this%phi_c , MPI_SUM, mesh%comm_parts)
        call XMPI_Allreduce(phi_d , this%phi_d , MPI_SUM, mesh%comm_parts)
        call XMPI_Allreduce(phi_ds, this%phi_ds, MPI_SUM, mesh%comm_parts)
        call XMPI_Allreduce(phi_d0, this%phi_d0, MPI_SUM, mesh%comm_parts)
        call XMPI_Allreduce(phi_s , this%phi_s , MPI_SUM, mesh%comm_parts)
        if (this % normalized) then
          this%phi_c  = this%phi_c  / volume
          this%phi_d  = this%phi_d  / volume
          this%phi_ds = this%phi_ds / volume
          this%phi_d0 = this%phi_d0 / volume
          this%phi_s  = this%phi_s  / volume
        end if
      end if

      ! finalization ...........................................................

      ! re-initialize parameters
      dx_min_loc = huge(ONE)
      dx_max_loc = 0
      v_max      = 0
      vv_max     = 0
      e_kin      = 0
      phi_c      = 0
      phi_d      = 0
      phi_ds     = 0
      phi_d0     = 0
      phi_s      = 0

      deallocate(sp, vp, mm, s0, s1, w)

      if (this % has_dissipation) then
        deallocate(bv_x, bv_u)
      end if

      !$omp end master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    end associate

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print header for incompressible flow characteristics

  subroutine PrintHeader(this, tag)
    class(INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    character(len=*), parameter :: fmt = '(1X,A6,4X)' !! short
  ! character(len=*), parameter :: fmt = '(2X,A6,5X)' !! long

    if (this % part == 0) then

      !$omp master
      write(*,'(A)',advance='NO') '#'
      write(*,fmt,advance='NO') '  t   '
      if (this % dt >= 0) then
        write(*,fmt,advance='NO') '  dt  '
      end if
      write(*,fmt,advance='NO') 'dx_min'
      write(*,fmt,advance='NO') 'dx_max'
      write(*,'(1X,A,2X)',advance='NO') 'po'
      write(*,fmt,advance='NO') 'v_max '
      write(*,fmt,advance='NO') 'e_kin '
      write(*,fmt,advance='NO') 'div_v '

      if (this % has_errors) then
        write(*,fmt,advance='NO') 'err_v '
        write(*,fmt,advance='NO') 'err_p '
      end if

      if (this % has_dissipation) then
        write(*,fmt,advance='NO') 'phi_c '
        write(*,fmt,advance='NO') 'phi_d '
        write(*,fmt,advance='NO') 'phi_ds'
        write(*,fmt,advance='NO') 'phi_d0'
        write(*,fmt,advance='NO') 'phi_s '
      end if

      if (present(tag)) then
        write(*,'(A)') tag
      else
        write(*,*)
      end if
      !$omp end master

    end if

  end subroutine PrintHeader

  !-----------------------------------------------------------------------------
  !> Print incompressible flow characteristics

  subroutine PrintValues(this, tag)
    class(INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    character(len=*), parameter :: fmt = '(ES10.3,1X)' !! short
  ! character(len=*), parameter :: fmt = '(ES12.5,1X)' !! long

    if (this % part == 0) then

      !$omp master
      write(*,fmt,advance='NO') this % t
      if (this % dt >= 0) then
        write(*,fmt,advance='NO') this % dt
      end if
      write(*,fmt,advance='NO') this % dx_min
      write(*,fmt,advance='NO') this % dx_max
      write(*,'(I4,1X)',advance='NO') this % po
      write(*,fmt,advance='NO') this % v_max
      write(*,fmt,advance='NO') this % e_kin
      write(*,fmt,advance='NO') this % div_v

      if (this % has_errors) then
        write(*,fmt,advance='NO') this % err_v
        write(*,fmt,advance='NO') this % err_p
      end if

      if (this % has_dissipation) then
        write(*,fmt,advance='NO') this % phi_c
        write(*,fmt,advance='NO') this % phi_d
        write(*,fmt,advance='NO') this % phi_ds
        write(*,fmt,advance='NO') this % phi_d0
        write(*,fmt,advance='NO') this % phi_s
      end if

      if (present(tag)) then
        write(*,'(X,A)') tag
      else
        write(*,*)
      end if
     !$omp end master

    end if

  end subroutine PrintValues

  !=============================================================================

end module INS__Flow_Characteristics__3D
