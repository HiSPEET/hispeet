!> summary:  Multilevel array reductions
!> author:   Joerg Stiller
!> date:     2024/10/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Array_Reductions__3D
  use Kind_Parameters
  use Array_Reductions
  use XMPI
  use ML__Mesh_Variable__3D
  implicit none
  private

  public :: ML_ScalarProduct_3D
  public :: ML_WeightedScalarProduct_3D

contains

  !-----------------------------------------------------------------------------
  !> Unweighted multilevel scalar product
  !>
  !> Computes the scalar product of `a` and `b` over active elements of all mesh
  !> levels. If `leaf` is present and true, the product is restricted to the
  !> leaf elements.
  !> The operands `a` and `b` are mesh variables with one ore more components
  !> per point. The result `r` is formed as the sum over the products of all
  !> components.
  !>
  !> Use `l_top` to specify a top level lower than `size(a%level)`

  real(RNP) function ML_ScalarProduct_3D(a, b, leaf, l_top) result(r)
    class(ML_MeshVariable_3D), intent(in) :: a
    class(ML_MeshVariable_3D), intent(in) :: b
    logical, optional, intent(in) :: leaf
    integer, optional, intent(in) :: l_top

    real(RNP), save :: r_glob, r_loc

    real(RNP) :: r_tmp
    logical   :: leaf_only
    integer   :: l_top_
    integer   :: i, l

    if (present(leaf)) then
      leaf_only = leaf
    else
      leaf_only = .false.
    end if

    if (present(l_top)) then
      l_top_ = min(l_top, size(a % level))
    else
      l_top_ = size(a % level)
    end if

    r_loc = 0

    do l = 1, l_top_
      associate( mesh => a % level(l) % mesh                  &
               , ne   => a % level(l) % mesh % n_elem_active  &
               , po   => a % level(l) % po                    &
               , nc   => a % level(l) % nc                    &
               , a_l  => a % level(l) % val                   &
               , b_l  => b % level(l) % val                   )

        if (l < l_top_ .and. leaf_only) then
          !$omp do reduction(+:r_loc) private(i)
          do i = 1, ne
            if (mesh % element(i) % IsLeaf()) then
              r_loc = r_loc + sum( a_l(0:po,0:po,0:po,i,1:nc) &
                                 * b_l(0:po,0:po,0:po,i,1:nc) )
            end if
          end do
        else
          r_tmp = 0
          do i = 1, nc
            r_tmp = r_tmp + ScalarProduct(a_l(:,:,:,:ne,i), b_l(:,:,:,:ne,i))
          end do
          !$omp master
          r_loc = r_loc + r_tmp
          !$omp end master
        end if

        if (l == l_top_) then
          !$omp master
          call XMPI_Allreduce(r_loc, r_glob, MPI_SUM, mesh%comm_world)
          r_loc = 0
          !$omp end master
          !$omp barrier
        end if

      end associate
    end do

    r = r_glob

  end function ML_ScalarProduct_3D

  !-----------------------------------------------------------------------------
  !> Weighted multilevel scalar product  --> separate module
  !>
  !> Computes the weighted scalar product of `a` and `b` over active elements of
  !> mesh levels. If `leaf` is present and true, the product is restricted to the
  !> leaf elements.
  !> The operands `a` and `b` are mesh variables with one ore more components
  !> per point. The weights `w` are associated with the mesh points and have
  !> either the same number of components ore only one.
  !>
  !> Use `l_top` to specify a top level lower than `size(a%level)`

  real(RNP) function ML_WeightedScalarProduct_3D(w, a, b, leaf, l_top) result(r)
    class(ML_MeshVariable_3D), intent(in) :: w
    class(ML_MeshVariable_3D), intent(in) :: a
    class(ML_MeshVariable_3D), intent(in) :: b
    logical, optional, intent(in) :: leaf
    integer, optional, intent(in) :: l_top

    real(RNP), save :: r_glob, r_loc

    real(RNP) :: r_tmp
    logical   :: leaf_only
    integer   :: l_top_
    integer   :: i, j, k, l, nw

    if (present(leaf)) then
      leaf_only = leaf
    else
      leaf_only = .false.
    end if

    if (present(l_top)) then
      l_top_ = min(l_top, size(a % level))
    else
      l_top_ = size(a % level)
    end if

    !$omp master
    r_loc = 0
    !$omp end master

    do l = 1, l_top_
      associate( mesh => a % level(l) % mesh                  &
               , ne   => a % level(l) % mesh % n_elem_active  &
               , po   => a % level(l) % po                    &
               , nc   => a % level(l) % nc                    &
               , a_l  => a % level(l) % val                   &
               , b_l  => b % level(l) % val                   &
               , w_l  => w % level(l) % val                   )

        nw = size(w_l, 5)

        if (l < l_top_ .and. leaf_only) then
          !$omp do collapse(2) reduction(+:r_loc) private(i,k)
          do k = 1, nc
          do i = 1, ne
            if (mesh % element(i) % IsLeaf()) then
              r_loc = r_loc + sum( w_l(0:po,0:po,0:po,i,min(k,nw)) &
                                 * a_l(0:po,0:po,0:po,i,k)         &
                                 * b_l(0:po,0:po,0:po,i,k)         )
            end if
          end do
          end do
        else
          r_tmp = 0
          do i = 1, nc
            j = min(i,nw)
            r_tmp = r_tmp + WeightedScalarProduct( w_l(:,:,:,:ne,j) &
                                                 , a_l(:,:,:,:ne,i) &
                                                 , b_l(:,:,:,:ne,i) )
          end do
          !$omp master
          r_loc = r_loc + r_tmp
          !$omp end master
        end if

      end associate
    end do

    !$omp master
    call XMPI_Allreduce(r_loc, r_glob, MPI_SUM, a%level(1)%mesh%comm_world)
    !$omp end master
    !$omp barrier

    r = r_glob

  end function ML_WeightedScalarProduct_3D

  !=============================================================================

end module ML__Array_Reductions__3D
