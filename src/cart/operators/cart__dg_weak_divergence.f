!> summary:  Weak divergence of a vector field
!> author:   Joerg Stiller
!> date:     2018/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__DG_Weak_Divergence

  use Kind_Parameters, only: RNP
  use Standard_Operators_1D
  use TPO__Div_3d_R
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Normal_Trace_Operator

  implicit none
  private

  public :: WeakDivergence

  interface WeakDivergence
    module procedure WeakDivergence_B
    module procedure WeakDivergence_O
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Weak divergence of a vector field on discontinuous elements - with BC

  subroutine WeakDivergence_B(mesh, Ms, Ds, u, bv_u, div_u)
    class(MeshPartition),    intent(in)  :: mesh           !< mesh partition
    real(RNP),               intent(in)  :: Ms(:)          !< std mass matrix
    real(RNP),               intent(in)  :: Ds(:,:)        !< std diff matrix
    real(RNP),               intent(in)  :: u(:,:,:,:,:)   !< vector field
    class(BoundaryVariable), intent(in)  :: bv_u(:)        !< BC
    real(RNP),               intent(out) :: div_u(:,:,:,:) !< weak divergence

    call WeakDivergence_X(mesh, size(Ms)-1, size(u,4), Ms, Ds, u, bv_u, div_u)

  end subroutine WeakDivergence_B

  !-----------------------------------------------------------------------------
  !> Weak divergence of a vector field on discontinuous elements - open (no BC)

  subroutine WeakDivergence_O(mesh, Ms, Ds, u, div_u)
    class(MeshPartition), intent(in)  :: mesh           !< mesh partition
    real(RNP),            intent(in)  :: Ms(:)          !< std mass matrix
    real(RNP),            intent(in)  :: Ds(:,:)        !< std diff matrix
    real(RNP),            intent(in)  :: u(:,:,:,:,:)   !< vector field
    real(RNP),            intent(out) :: div_u(:,:,:,:) !< weak divergence

    call WeakDivergence_X(mesh, size(Ms)-1, size(u,4), Ms, Ds, u, div_u=div_u)

  end subroutine WeakDivergence_O

  !-----------------------------------------------------------------------------
  !> Weak divergence of a vector field on discontinuous elements - eXplicit

  subroutine WeakDivergence_X(mesh, po, ne, Ms, Ds, u, bv_u, div_u)

    ! arguments ................................................................

    !> mesh partition
    class(MeshPartition), intent(in) :: mesh
    !> polynomial order
    integer,   intent(in) :: po
    !> number of elements
    integer,   intent(in) :: ne
    !> std mass matrix
    real(RNP), intent(in) :: Ms(0:po)
    !> std diff matrix
    real(RNP), intent(in) :: Ds(0:po, 0:po)
    !> vector field
    real(RNP), intent(in) :: u(0:po, 0:po, 0:po, ne, 3)
    !> boundary conditions for u
    class(BoundaryVariable), optional, intent(in) :: bv_u(mesh%n_boundary)
    !> weak divergence of u
    real(RNP), intent(out) :: div_u(0:po, 0:po, 0:po, ne)

    ! internal variables .......................................................

    type(NormalTraceOperator), asynchronous, allocatable, save :: trace_op
    real(RNP), allocatable, save :: tr_un(:,:,:,:)
    real(RNP) :: g, qf(0:po,0:po)

    integer :: tag = 1000
    integer :: e, f, f1, f2

    ! initialization ...........................................................

    !$omp single
    allocate(trace_op)
    allocate(tr_un(0:po, 0:po, 2, mesh%nf))
    !$omp end single

    ! extract and start transferring traces ....................................

    if (present(bv_u)) then
      call trace_op % GetTrace_Start(mesh, u, bv_u, tr_un, tag)
    else
      call trace_op % GetTrace_Start(mesh, u, tr_un, tag)
    end if

    ! compute transposed divergence ............................................

    call TPO_Div_R(Ds, mesh%dx, u, div_u)

    ! complete transfer ........................................................

    call trace_op % GetTrace_Finish(mesh, tr_un)

    associate(dx => mesh%dx, face => mesh%face)

      ! x1-flux contribution ...................................................

      g  = -1 / (dx(1) * Ms(0))   ! metric and averaging factor
      f1 = 1                      ! first face
      f2 = mesh % nf1             ! last

      !$omp do
      do f = f1, f2

        qf = g * (tr_un(:,:,1,f) + tr_un(:,:,2,f))

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne) then
          div_u(po,:,:,e) = div_u(po,:,:,e) + qf
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne) then
          div_u(0,:,:,e) = div_u(0,:,:,e) + qf
        end if

      end do

      ! x2-flux contribution ...................................................

      g  = -1 / (dx(2) * Ms(0))   ! metric and averaging factor
      f1 = f2 + 1                 ! first face
      f2 = f2 + mesh % nf2        ! last

      !$omp do
      do f = f1, f2

        qf = g * (tr_un(:,:,1,f) + tr_un(:,:,2,f))

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne) then
          div_u(:,po,:,e) = div_u(:,po,:,e) + qf
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne) then
          div_u(:,0,:,e) = div_u(:,0,:,e) + qf
        end if

      end do

      ! x3-flux contribution ...................................................

      g  = -1 / (dx(3) * Ms(0))   ! metric and averaging factor
      f1 = f2 + 1                 ! first face
      f2 = f2 + mesh % nf3        ! last

      !$omp do
      do f = f1, f2

        qf = g * (tr_un(:,:,1,f) + tr_un(:,:,2,f))

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne) then
          div_u(:,:,po,e) = div_u(:,:,po,e) + qf
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne) then
          div_u(:,:,0,e) = div_u(:,:,0,e) + qf
        end if

      end do

    end associate

    ! finalization .............................................................

    !$omp single
    deallocate(tr_un)
    deallocate(trace_op)
    !$omp end single

  end subroutine WeakDivergence_X

  !=============================================================================

end module CART__DG_Weak_Divergence
