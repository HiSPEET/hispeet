!> summary:  3D mesh metrics
!> author:   Joerg Stiller
!> date:     2021/05/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - more efficient (separate) treatment of regular meshes
!===============================================================================

module Mesh_Metrics__3D
  use Kind_Parameters, only: RNP
  use Standard_Element_Operators__1D
  use Mesh__3D
  implicit none
  private

  public :: MeshMetrics_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh metrics
  !>
  !> Provides the mesh points and metric coefficients for all elements of the
  !> given partition.
  !>
  !> ### Components referring to points 'i,j,k' in element `e`
  !>
  !>   -  `x  (i,j,k,e,:)  `    : Cartesian coordinates
  !>   -  `Jm (i,j,k,e,:,:)`    : Jacobian matrix J = ∂x/∂ξ
  !>   -  `Ji (i,j,k,e,:,:)`    : inverse Jacobian matrix J⁻¹ = ∂ξ/∂x
  !>   -  `Jd (i,j,k,e)    `    : determinant of the Jacobian matrix ‖J‖
  !>   -  `G  (i,j,k,e,:)  `    : Laplacian metrics  ‖J‖J⁻¹J⁻ᵀ
  !>
  !> As G is symmetric, only its upper triangular part is stored, i.e:
  !>
  !>       G(... , 1) = G₁₁
  !>       G(... , 2) = G₁₂
  !>       G(... , 3) = G₁₃
  !>       G(... , 4) = G₂₂
  !>       G(... , 5) = G₂₃
  !>       G(... , 6) = G₃₃
  !>
  !> ### Components referring to points 'l,m' on face `f` of element `e`
  !>
  !>   -  `a      (l,m,f,e)  `  : area coefficient a =|∂x/∂τ₁ × ∂x/∂τ₂|
  !>   -  `n      (l,m,f,e,:)`  : unit normal vector n = (∂x/∂τ₁ × ∂x/∂τ₂)/a
  !>   -  `Ji_n_a (l,m,f,e,:)`  : J⁻¹⋅n

  type MeshMetrics_3D
    real(RNP), allocatable :: x    (:,:,:,:,:)   !< mesh points
    real(RNP), allocatable :: Jm   (:,:,:,:,:,:) !< Jacobian matrix
    real(RNP), allocatable :: Ji   (:,:,:,:,:,:) !< Jacobian matrix inverse
    real(RNP), allocatable :: Jd   (:,:,:,:)     !< Jacobian matrix determinant
    real(RNP), allocatable :: G    (:,:,:,:,:)   !< Laplacian metrics
    real(RNP), allocatable :: a    (:,:,:,:)     !< area coefficient
    real(RNP), allocatable :: n    (:,:,:,:,:)   !< unit normal vector
    real(RNP), allocatable :: Ji_n (:,:,:,:,:)   !< J⁻¹⋅n
  contains
    procedure :: Init_MeshMetrics_3D
  end type MeshMetrics_3D

  ! constructor
  interface MeshMetrics_3D
    module procedure New_MeshMetrics_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> MeshMetrics_3D constructor

  function New_MeshMetrics_3D(mesh, std_op) result(this)
    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    class(StandardElementOperators_1D), intent(in) :: std_op !< standard ops

    type(MeshMetrics_3D) :: this

    call Init_MeshMetrics_3D(this, mesh, std_op)

  end function New_MeshMetrics_3D

  !-----------------------------------------------------------------------------
  !> Initialization of mesh metrics

  subroutine Init_MeshMetrics_3D(this, mesh, std_op)
    class(MeshMetrics_3D), intent(inout) :: this
    class(Mesh_3D),        intent(in)    :: mesh !< mesh partition
    class(StandardElementOperators_1D), intent(in) :: std_op !< standard ops

    ! local variables ..........................................................

    real(RNP), allocatable :: Dt(:,:), grad_x(:,:,:,:,:)
    real(RNP) :: a, c, Jd, Jm(3,3), Ji(3,3), n(3)
    integer   :: e, f, i, j, k, p, q, np, ne

    associate(po => std_op % po)

      ! basic initialization ...................................................

      np = po + 1
      ne = mesh % n_elem

      call mesh % GetPoints(po, std_op % basis, this % x)

      allocate( this % Jm (0:po,0:po,0:po,ne,3,3) )
      allocate( this % Ji (0:po,0:po,0:po,ne,3,3) )
      allocate( this % Jd (0:po,0:po,0:po,ne)     )
      allocate( this % G  (0:po,0:po,0:po,ne,6)   )

      allocate( this % a    (0:po,0:po,6,ne)   )
      allocate( this % n    (0:po,0:po,6,ne,3) )
      allocate( this % Ji_n (0:po,0:po,6,ne,3) )

      ! auxiliary data
      allocate(Dt(0:po,0:po), grad_x(0:po,0:po,0:po,3,3))

      Dt = transpose(std_op % D)

      !$omp do
      do e = 1, ne

        ! grad_x = dx/dξ .......................................................

        call TPO_Grad_X(np, Dt, this % x(:,:,:,e,1), grad_x(:,:,:,:,1))
        call TPO_Grad_X(np, Dt, this % x(:,:,:,e,2), grad_x(:,:,:,:,2))
        call TPO_Grad_X(np, Dt, this % x(:,:,:,e,3), grad_x(:,:,:,:,3))

        ! volume metrics .......................................................

        do k = 0, po
        do j = 0, po
        do i = 0, po

          ! Jacobian matrix
          do q = 1, 3
          do p = 1, 3
            Jm(p,q) = grad_x(i,j,k,q,p)
          end do
          end do

          ! Jacobian determinant
          Jd = Jm(1,1) * (Jm(2,2)*Jm(3,3) - Jm(2,3)*Jm(3,2))  &
             + Jm(1,2) * (Jm(2,3)*Jm(3,1) - Jm(2,1)*Jm(3,3))  &
             + Jm(1,3) * (Jm(2,1)*Jm(3,2) - Jm(2,2)*Jm(3,1))

          ! inverse Jacobian matrix
          c = 1 / max(Jd, epsilon(c))

          Ji(1,1) = c * (Jm(2,2)*Jm(3,3) - Jm(3,2)*Jm(2,3))
          Ji(2,1) = c * (Jm(2,3)*Jm(3,1) - Jm(3,3)*Jm(2,1))
          Ji(3,1) = c * (Jm(2,1)*Jm(3,2) - Jm(3,1)*Jm(2,2))

          Ji(1,2) = c * (Jm(3,2)*Jm(1,3) - Jm(1,2)*Jm(3,3))
          Ji(2,2) = c * (Jm(3,3)*Jm(1,1) - Jm(1,3)*Jm(3,1))
          Ji(3,2) = c * (Jm(3,1)*Jm(1,2) - Jm(1,1)*Jm(3,2))

          Ji(1,3) = c * (Jm(1,2)*Jm(2,3) - Jm(2,2)*Jm(1,3))
          Ji(2,3) = c * (Jm(1,3)*Jm(2,1) - Jm(2,3)*Jm(1,1))
          Ji(3,3) = c * (Jm(1,1)*Jm(2,2) - Jm(2,1)*Jm(1,2))

          ! assignment
          this % Jd(i,j,k,e) = Jd
          do q = 1, 3
          do p = 1, 3
            this % Jm(i,j,k,e,p,q) = Jm(p,q)
            this % Ji(i,j,k,e,p,q) = Ji(p,q)
          end do
          end do
          this % G (i,j,k,e,1) = Jd * sum(Ji(1,:) * Ji(1,:))
          this % G (i,j,k,e,2) = Jd * sum(Ji(1,:) * Ji(2,:))
          this % G (i,j,k,e,3) = Jd * sum(Ji(1,:) * Ji(3,:))
          this % G (i,j,k,e,4) = Jd * sum(Ji(2,:) * Ji(2,:))
          this % G (i,j,k,e,5) = Jd * sum(Ji(2,:) * Ji(3,:))
          this % G (i,j,k,e,6) = Jd * sum(Ji(3,:) * Ji(3,:))

        end do
        end do
        end do

        ! face 1: ξ₁ = -1 ......................................................

        f = 1
        i = 0
        do k = 0, po
        do j = 0, po

          n = - this % Jd(i,j,k,e) * this % Ji(i,j,k,e,1,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(j,k,f,e)   = a
          this % n(j,k,f,e,:) = n

          this % Ji_n(j,k,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

        ! face 2: ξ₁ =  1 ......................................................

        f = 2
        i = po
        do k = 0, po
        do j = 0, po

          n = + this % Jd(i,j,k,e) * this % Ji(i,j,k,e,1,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(j,k,f,e)   = a
          this % n(j,k,f,e,:) = n

          this % Ji_n(j,k,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

        ! face 3: ξ₂ = -1 ......................................................

        f = 3
        j = 0
        do k = 0, po
        do i = 0, po

          n = - this % Jd(i,j,k,e) * this % Ji(i,j,k,e,2,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(i,k,f,e)   = a
          this % n(i,k,f,e,:) = n

          this % Ji_n(i,k,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

        ! face 4: ξ₂ =  1 ......................................................

        f = 4
        j = po
        do k = 0, po
        do i = 0, po

          n = + this % Jd(i,j,k,e) * this % Ji(i,j,k,e,2,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(i,k,f,e)   = a
          this % n(i,k,f,e,:) = n

          this % Ji_n(i,k,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

        ! face 5: ξ₃ = -1 ......................................................

        f = 5
        k = 0
        do j = 0, po
        do i = 0, po

          n = - this % Jd(i,j,k,e) * this % Ji(i,j,k,e,3,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(i,j,f,e)   = a
          this % n(i,j,f,e,:) = n

          this % Ji_n(i,j,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

        ! face 6: ξ₃ =  1 .......................................................

        f = 6
        k = po
        do j = 0, po
        do i = 0, po

          n = + this % Jd(i,j,k,e) * this % Ji(i,j,k,e,3,:)
          a = sqrt(sum(n**2))
          n = n / a

          this % a(i,j,f,e)   = a
          this % n(i,j,f,e,:) = n

          this % Ji_n(i,j,f,e,:) = this % Ji(i,j,k,e,:,1) * n(1) &
                                 + this % Ji(i,j,k,e,:,2) * n(2) &
                                 + this % Ji(i,j,k,e,:,3) * n(3)
        end do
        end do

      end do

      ! finalization ...........................................................

      deallocate(grad_x)

    end associate

  end subroutine Init_MeshMetrics_3D

  !-----------------------------------------------------------------------------
  !> Standard element gradient

  pure subroutine TPO_Grad_X(np, Dt, x, g)
    integer,   intent(in)  :: np            !< num points per direction
    real(RNP), intent(in)  :: Dt(np,np)     !< transposed standard diff matrix
    real(RNP), intent(in)  :: x(np,np,np)   !< scalar element variable
    real(RNP), intent(out) :: g(np,np,np,3) !< gradient of x

    integer :: i, j, k, p

    ! v1 = du/dx1, v2 = du/dx2 .................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      g(i,j,k,1) = 0
      g(i,j,k,2) = 0
      do p = 1, np
        g(i,j,k,1) = g(i,j,k,1) + Dt(p,i) * x(p,j,k)
        g(i,j,k,2) = g(i,j,k,2) + Dt(p,j) * x(i,p,k)
      end do
    end do
    end do
    end do

    ! v3 = du/dx3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      g(i,j,k,3) = 0
      do p = 1, np
        g(i,j,k,3) =  g(i,j,k,3) + Dt(p,k) * x(i,j,p)
      end do
    end do
    end do
    end do

  end subroutine TPO_Grad_X

  !=============================================================================

end module Mesh_Metrics__3D
