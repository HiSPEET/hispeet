!> summary:  3D spacetime solution projection from fine to coarse level
!> author:   Joerg Stiller
!> date:     2024/12/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module FC__Spacetime_Projection__3D
  use Constants, only: ZERO
  use HP__Coarsening_Operator__1D
  use Mesh_Variable__3D
  use Spacetime_Variable__3D
  use Child_To_Parent_Projection__3D
  use ML__Spacetime_Operators__3D

  implicit none
  private

  public :: FC_SpacetimeProjection_3D

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse space-time projection of solution-like variables
  !>
  !> Use `method` to select either
  !>   - `'I'`  embedded interpolation or
  !>   - `'P'`  L² projection

  subroutine FC_SpacetimeProjection_3D(ml_op, method, l_f, v_f, v_c)
    class(ML_SpacetimeOperators_3D), intent(in)    :: ml_op
    character,                       intent(in)    :: method !< 'I' or 'P'
    integer,                         intent(in)    :: l_f    !< fine level ID
    class(SpacetimeVariable_3D),     intent(in)    :: v_f    !< fine variable
    class(SpacetimeVariable_3D),     intent(inout) :: v_c    !< coarse variable

    type(SpacetimeVariable_3D), allocatable, save :: v_i
    type(HP_CoarseningOperator_1D), allocatable, save :: pop_x, pop_t
    integer :: l_c
    integer :: px_c, pt_c, nt_c
    integer :: pt_f, nt_f
    integer :: c, e, k, m, n, n1, n2, nc

    ! initialization ...........................................................

    l_c = l_f - 1

    ! coarse mesh dimensions
    px_c = ml_op % sem(l_c) % std_op % po
    pt_c = v_c % po_t
    nt_c = v_c % ne_t

    ! fine mesh dimensions
    pt_f = v_f % po_t
    nt_f = v_f % ne_t

    ! number of components
    nc = v_f % nc

    !$omp master
    select case(method)
    case('I')
      ! use embedded interpolation
      pop_x = ml_op % iop_fc_x(l_f)
      pop_t = ml_op % iop_fc_t(l_f)
    case('P')
      ! use L² projection
      pop_x = ml_op % pop_fc_x(l_f)
      pop_t = ml_op % pop_fc_t(l_f)
    end select
    allocate(v_i)
    v_i % po_t = pt_f
    v_i % ne_t = nt_f
    allocate(v_i % var(0:pt_f, 1:nt_f))
    do n = 1, nt_f
    do m = 0, pt_f
      v_i % var(m,n) = MeshVariable_3D(ml_op % sem(l_c) % mesh, px_c, nc)
    end do
    end do
    !$omp end master

    ! spatial projection .......................................................

    do n = 1, nt_f
    do m = 0, pt_f
      call ChildToParentProjection_3D( ml_op % sem(l_f) % mesh & ! child     c
                                     , ml_op % sem(l_c) % mesh & ! parent    p
                                     , pop_x                   & ! iop: c -> p
                                     , v_f % var(m,n) % val    & ! v_c
                                     , v_i % var(m,n) % val    ) ! v_p
    end do
    end do

    ! temporal projection ......................................................

    select case(pop_t % mode)

    case(0)

      ! identity: pt_c = pt_f, nt_c = nt_f
      do n = 1, nt_c
      do m = 0, pt_c
        associate( mesh => ml_op % sem(l_c) % mesh)
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, mesh % n_elem_active
            if (mesh % element(e) % adaptation % refinement < 1000) cycle
            v_c % var(m,n) % val(:,:,:,e,c) = v_i % var(m,n) % val(:,:,:,e,c)
          end do
          end do
        end associate
      end do
      end do

    case(1)

      ! p-coarsening: nt_c = nt_f
      do n = 1, nt_c
      do m = 0, pt_c
        associate( mesh => ml_op % sem(l_c) % mesh)
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, mesh % n_elem_active
            if (mesh % element(e) % adaptation % refinement < 1000) cycle
            v_c % var(m,n) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_f
              v_c % var(m,n) % val(:,:,:,e,c)                          &
                  = v_c % var(m,n) % val(:,:,:,e,c)                    &
                  + v_i % var(k,n) % val(:,:,:,e,c) * pop_t % A(m,k,1)
            end do
          end do
          end do
        end associate
      end do
      end do

    case(2)

      ! hp-coarsening: 2 * nt_c = nt_f
      ! use combined smoothing+projection operator pop_t % C
      do n = 1, nt_c
      do m = 0, pt_c
        n1 = 2 * n - 1
        n2 = 2 * n
        associate( mesh => ml_op % sem(l_c) % mesh)
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, mesh % n_elem_active
            if (mesh % element(e) % adaptation % refinement < 1000) cycle
            v_c % var(m,n) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_f
              v_c % var(m,n) % val(:,:,:,e,c)                            &
                  = v_c % var(m,n ) % val(:,:,:,e,c)                     &
                  + v_i % var(k,n1) % val(:,:,:,e,c) * pop_t % C(m,k,1)  &
                  + v_i % var(k,n2) % val(:,:,:,e,c) * pop_t % C(m,k,2)
            end do
          end do
          end do
        end associate
      end do
      end do

    end select

    ! finalization .............................................................

    !$omp master
    deallocate(pop_x, pop_t, v_i)
    !$omp end master

  end subroutine FC_SpacetimeProjection_3D

  !=============================================================================

end module FC__Spacetime_Projection__3D
