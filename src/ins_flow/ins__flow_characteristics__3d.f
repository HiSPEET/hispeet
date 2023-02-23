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
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_FlowCharacteristics_3D

  !-----------------------------------------------------------------------------
  !> Incompressible flow statistics

  type INS_FlowCharacteristics_3D

    real(RNP) :: t     = -1 !< problem time
    real(RNP) :: err_v = -1 !< L² velocity error
    real(RNP) :: err_p = -1 !< L² velocity error
    real(RNP) :: div_v = -1 !< L² velocity divergence
    real(RNP) :: v_max = -1 !< maximum velocity
    real(RNP) :: e_kin = -1 !< kinetic energy per unit mass

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

  subroutine Evaluate(this, problem, ins_op, t, u, volume)
    class(INS_FlowCharacteristics_3D), intent(inout) :: this
    class(INS_Problem_3D),  intent(in) :: problem
    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP),              intent(in) :: t
    real(RNP), contiguous,  intent(in) :: u(:,:,:,:,:)
    real(RNP), optional,    intent(in) :: volume

    real(RNP), allocatable, save :: vp(:,:,:,:,:), w(:,:,:,:,:)
    real(RNP), save :: v_max, vv_max = 0
    real(RNP), save :: e_kin, vv_sum = 0
    real(RNP) :: div_v, err_v, err_p, vv

    integer :: e, i, j, k, ne, np

    associate( mesh => ins_op % mesh                &
             , eop  => ins_op % eop_v               &
             , x    => ins_op % sem_v % metrics % x )

      ! initialization .........................................................

      if (mesh % part < 0) then
        !$omp barrier
        return
      end if

      np = size(u,1)
      ne = size(u,4)

      !$omp master

      this % t          = t
      this % part       = mesh % part
      this % has_errors = problem % HasExactSolution()

      allocate(vp(np,np,6,ne,3))
      allocate(w, mold = u)

      !$omp end master
      !$omp barrier

      ! errors .................................................................

      if (this % has_errors) then
        associate(q => w(:,:,:,:,4))

          call problem % GetExactSolution(x, t, w)        !
          call MergeArrays(-ONE, w, ONE, u, multi=.true.) ! w  = u  - u_ex
          call CalibrateArray(q, comm = mesh%comm_parts)  ! subtract avg p-error

          ! pressure
          !$omp do
          do e = 1, ne
            q(:,:,:,e) = q(:,:,:,e)**2
          end do
          call GetVolumeIntegral(ins_op%sem_v, q, err_p)
          err_p = sqrt(err_p)

          ! velocity
          !$omp do
          do e = 1, ne
            q(:,:,:,e) = w(:,:,:,e,1)**2 + w(:,:,:,e,2)**2 + w(:,:,:,e,3)**2
          end do
          call GetVolumeIntegral(ins_op%sem_v, q, err_v)
          err_v = sqrt(err_v)

        end associate
      end if

      ! divergence .............................................................

      associate(v => u(:,:,:,:,1:3), q => w(:,:,:,:,4))
        call GetOuterTraces_3D(mesh, v, vp)         ! vp = v⁺
        call TPO_Div(eop, ins_op%sem_v, v, vp, q)   ! q = ∇⋅v
        !$omp do
        do e = 1, ne
          q(:,:,:,e) = q(:,:,:,e) ** 2
        end do
        call GetVolumeIntegral(ins_op%sem_v, q, div_v)
        div_v = sqrt(div_v)
      end associate

      ! maximum velocity and energy ............................................

      !$omp do reduction(max:vv_max, sum:vv_sum)
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np
          vv = u(i,j,k,e,1)**2 + u(i,j,k,e,2)**2 + u(i,j,k,e,3)**2
          vv_max = max(vv_max, vv)
          vv_sum = vv_sum + vv
        end do
        end do
        end do
      end do

      !$omp master
      call XMPI_Allreduce(sqrt(vv_max), v_max, MPI_MAX, mesh%comm_parts)
      call XMPI_Allreduce(HALF*vv_sum , e_kin, MPI_SUM, mesh%comm_parts)
      !$omp end master

      ! finalization ...........................................................

      !$omp master

      if (present(volume)) then
        this % err_v = err_v / volume
        this % err_p = err_p / volume
        this % div_v = div_v / volume
        this % e_kin = e_kin / volume
      else
        this % err_v  = err_v
        this % err_p  = err_p
        this % div_v  = div_v
      end if
      this % v_max = v_max

      vv_max = 0
      vv_sum = 0

      deallocate(vp, w)

      !$omp end master

    end associate

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print header for incompressible flow characteristics

  subroutine PrintHeader(this)
    class(INS_FlowCharacteristics_3D), intent(in) :: this

    if (this % part == 0) then

      !$omp master
      write(*,'(A,2X)'   ,advance='NO') '#'
      write(*,'(2X,A,8X)',advance='NO') 't'
      write(*,'(2X,A,6X)',advance='NO') 'v_max'
      write(*,'(2X,A,6X)',advance='NO') 'e_kin'
      write(*,'(2X,A,6X)',advance='NO') 'div_v'

      if (this % has_errors) then
        write(*,'(2X,A,6X)',advance='NO') 'err_v'
        write(*,'(2X,A,6X)',advance='NO') 'err_p'
      end if

      write(*,'(A)') '#char#'
      !$omp end master

    end if

  end subroutine PrintHeader

  !-----------------------------------------------------------------------------
  !> Print incompressible flow characteristics

  subroutine PrintValues(this)
    class(INS_FlowCharacteristics_3D), intent(in) :: this

    if (this % part == 0) then

      !$omp master
      write(*,'(ES12.5,1X)',advance='NO') this % t
      write(*,'(ES12.5,1X)',advance='NO') this % v_max
      write(*,'(ES12.5,1X)',advance='NO') this % e_kin
      write(*,'(ES12.5,1X)',advance='NO') this % div_v

      if (this % has_errors) then
        write(*,'(ES12.5,1X)',advance='NO') this % err_v
        write(*,'(ES12.5,1X)',advance='NO') this % err_p
      end if

      write(*,'(1X,A)') '#char#'
      !$omp end master

    end if

  end subroutine PrintValues

  !=============================================================================

end module INS__Flow_Characteristics__3D
