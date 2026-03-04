!> summary:  Evaluation of flow characteristics
!> author:   Joerg Stiller
!> date:     2023/02/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Flow_Characteristics__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use Array_Reductions
  use XMPI

  use Trace_Operators__3D
  use TPO__Div__3D
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

    integer, private :: part = -1
    logical, private :: has_errors = .false.

  contains
    procedure :: Evaluate
    procedure :: PrintHeader
    procedure :: PrintValues
  end type INS_FlowCharacteristics_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of incompressible flow characteristics

  subroutine Evaluate(this, ins_op, t, u, dt, volume, leaf)
    class(INS_FlowCharacteristics_3D), intent(inout) :: this
    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP),              intent(in) :: t
    real(RNP), contiguous,  intent(in) :: u(:,:,:,:,:)
    real(RNP), optional,    intent(in) :: dt
    real(RNP), optional,    intent(in) :: volume
    logical,   optional,    intent(in) :: leaf !< constrain to leaf elements [F]

    real(RNP), allocatable, save :: vp(:,:,:,:,:), w(:,:,:,:,:)
    real(RNP), save :: v_max, vv_max
    real(RNP), save :: e_kin
    real(RNP), save :: dx_min_loc, dx_max_loc
    real(RNP) :: div_v, err_v, err_p, vv
    real(RNP) :: dx(3)

    logical :: complete
    integer :: e, i, j, k, na, nc, ne, np

    if (present(leaf)) then
      complete = .not. leaf
    else
      complete = .true.
    end if

    associate( problem => ins_op % problem             &
             , mesh    => ins_op % mesh                &
             , eop     => ins_op % eop_u               &
             , x       => ins_op % sem_u % metrics % x )

      ! initialization .........................................................

      if (mesh % part < 0) then
        !$omp barrier
        return
      end if

      np = size(u,1)
      ne = size(u,4)
      nc = size(u,5)
      na = mesh % n_elem_active

      !$omp master

      this % t          = t
      this % part       = mesh % part
      this % has_errors = problem % HasExactSolution()

      if (present(dt)) then
        this % dt = dt
      else
        this % dt = -1
      end if

      this % po     = ins_op % eop_u % po
      this % err_v  = -1
      this % err_p  = -1

      dx_min_loc = huge(ONE)
      dx_max_loc = 0
      vv_max = 0

      allocate(vp(np,np,6,ne,3))
      allocate(w(np,np,np,ne,nc), source = ZERO)

      !$omp end master
      !$omp barrier

      ! mesh spacing ...........................................................

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

      ! errors .................................................................

      if (this % has_errors) then
        associate(q => w(:,:,:,:,4))

          call problem % GetExactSolution(x, t, w)         !
          call MergeArrays(-ONE, w, ONE, u, multi=.true.)  ! w  = u  - u_ex
          if (complete .and. ne == na) then
            ! remove mean pressure error
            call CalibrateArray(q, comm = mesh%comm_parts)
          end if

          ! pressure
          !$omp do
          do e = 1, na
            q(:,:,:,e) = q(:,:,:,e)**2
          end do
          call GetVolumeIntegral(ins_op%sem_u, q, err_p, leaf)
          err_p = sqrt(err_p)

          ! velocity
          !$omp do
          do e = 1, na
            q(:,:,:,e) = w(:,:,:,e,1)**2 + w(:,:,:,e,2)**2 + w(:,:,:,e,3)**2
          end do
          call GetVolumeIntegral(ins_op%sem_u, q, err_v, leaf)
          err_v = sqrt(err_v)

        end associate
      end if

      ! divergence .............................................................

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

      ! maximum velocity and energy ............................................

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

      ! finalization ...........................................................

      !$omp master

      if (present(volume)) then
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

      vv_max = 0

      deallocate(vp, w)

      !$omp end master

    end associate

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print header for incompressible flow characteristics

  subroutine PrintHeader(this, tag)
    class(INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    if (this % part == 0) then

      !$omp master
      write(*,'(A,2X)'   ,advance='NO') '#'
      write(*,'(3X,A,8X)',advance='NO') 't'
      if (this % dt >= 0) then
        write(*,'(3X,A,7X)',advance='NO') 'dt'
      end if
      write(*,'(2X,A,5X)',advance='NO') 'dx_min'
      write(*,'(2X,A,5X)',advance='NO') 'dx_max'
      write(*,'(1X,A,2X)',advance='NO') 'po'
      write(*,'(2X,A,6X)',advance='NO') 'v_max'
      write(*,'(2X,A,6X)',advance='NO') 'e_kin'
      write(*,'(2X,A,6X)',advance='NO') 'div_v'

      if (this % has_errors) then
        write(*,'(2X,A,6X)',advance='NO') 'err_v'
        write(*,'(2X,A,6X)',advance='NO') 'err_p'
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

    if (this % part == 0) then

      !$omp master
      write(*,'(ES12.5,1X)',advance='NO') this % t
      if (this % dt >= 0) then
        write(*,'(ES12.5,1X)',advance='NO') this % dt
      end if
      write(*,'(ES12.5,1X)',advance='NO') this % dx_min
      write(*,'(ES12.5,1X)',advance='NO') this % dx_max
      write(*,'(I4    ,1X)',advance='NO') this % po
      write(*,'(ES12.5,1X)',advance='NO') this % v_max
      write(*,'(ES12.5,1X)',advance='NO') this % e_kin
      write(*,'(ES12.5,1X)',advance='NO') this % div_v

      if (this % has_errors) then
        write(*,'(ES12.5,1X)',advance='NO') this % err_v
        write(*,'(ES12.5,1X)',advance='NO') this % err_p
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
