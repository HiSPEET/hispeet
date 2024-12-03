!> summary:  3D spacetime residual restriction from fine to coarse level
!> author:   Joerg Stiller
!> date:     2024/12/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module FC__Spacetime_Restriction__3D
  use Constants, only: ZERO
  use Array_Assignments
  use Mesh_Variable__3D
  use Spacetime_Variable__3D
  use Child_To_Parent_Restriction__3D
  use ML__Spacetime_Operators__3D

  implicit none
  private

  public :: FC_SpacetimeRestriction_3D

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse space-time restriction of residual-like variables

  subroutine FC_SpacetimeRestriction_3D(ml_op, l_f, v_f, v_c)
    class(ML_SpacetimeOperators_3D), intent(in)    :: ml_op
    integer,                         intent(in)    :: l_f !< fine level ID
    class(SpacetimeVariable_3D),     intent(in)    :: v_f !< fine variable
    class(SpacetimeVariable_3D),     intent(inout) :: v_c !< coarse variable

    type(SpacetimeVariable_3D), allocatable, save :: v_i
    integer :: l_c
    integer :: px_c, nx_c, pt_c, nt_c
    integer :: px_f, nx_f, pt_f, nt_f
    integer :: c, e, k, m, n, n1, n2, nc

    ! initialization ...........................................................

    l_c = l_f - 1

    ! coarse mesh dimensions
    px_c = ml_op % sem(l_c) % std_op % po
    nx_c = ml_op % sem(l_c) % mesh % n_elem
    pt_c = v_c % po_t
    nt_c = v_c % ne_t

    ! fine mesh dimensions
    px_f = ml_op % sem(l_f) % std_op % po
    nx_f = ml_op % sem(l_f) % mesh % n_elem
    pt_f = v_f % po_t
    nt_f = v_f % ne_t

    ! number of components
    nc = v_f % var(0,1) % nc

    !$omp master
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

    ! spatial restriction ......................................................

    do n = 1, nt_f
    do m = 0, pt_f
      call ChildToParentRestriction_3D( ml_op % sem(l_f) % mesh & ! child     c
                                      , ml_op % sem(l_c) % mesh & ! parent    p
                                      , ml_op % iop_cf_x(l_c)   & ! iop: c -> p
                                      , v_f % var(m,n) % val    & ! v_c
                                      , v_i % var(m,n) % val    ) ! v_p
    end do
    end do

    ! temporal interpolation ...................................................

    associate(A => ml_op % iop_cf_t(l_c) % A)

      select case(size(A,3))

      case(0)

        ! identity: pt_c = pt_f, nt_c = nt_f
        do n = 1, nt_c
        do m = 0, pt_c
          call SetArray( v_f % var(m,n) % val &
                       , v_i % var(m,n) % val &
                       , multi = .true.       )
        end do
        end do

      case(1)

        ! p-coarsening: nt_c = nt_f
        do n = 1, nt_c
        do m = 0, pt_c
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, nx_c
            v_c % var(m,n) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_f
              v_c % var(m,n) % val(:,:,:,e,c)                   &
                  = v_c % var(m,n) % val(:,:,:,e,c)             &
                  + v_i % var(k,n) % val(:,:,:,e,c) * A(k,m,1)
            end do
          end do
          end do
        end do
        end do

      case(2)

        ! hp-coarsening: 2 * nt_c = nt_f
        do n = 1, nt_c
        do m = 0, pt_c
          n1 = 2 * n - 1
          n2 = 2 * n
          !$omp do collapse(2)
          do c = 1, nc
          do e = 1, nx_c
            v_c % var(m,n) % val(:,:,:,e,c) = ZERO
            do k = 0, pt_f
              v_c % var(m,n) % val(:,:,:,e,c)                    &
                  = v_c % var(m,n ) % val(:,:,:,e,c)             &
                  + v_i % var(k,n1) % val(:,:,:,e,c) * A(k,m,1)  &
                  + v_i % var(k,n2) % val(:,:,:,e,c) * A(k,m,2)
            end do
          end do
          end do
        end do
        end do

      end select

    end associate

    ! finalization .............................................................

    !$omp master
    deallocate(v_i)
    !$omp end master

  end subroutine FC_SpacetimeRestriction_3D

  !=============================================================================

end module FC__Spacetime_Restriction__3D
