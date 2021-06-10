!> summary:  3D mesh metrics
!> author:   Joerg Stiller
!> date:     2021/05/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   *  replace `TPO_Grad_S` by optimized external procedure
!>   *  OpenMP parallelization
!===============================================================================

module Mesh_3d__Metrics
  use Kind_Parameters, only: RNP
  use Standard_Operators__1D
  use Mesh_3d__Partition
  implicit none
  private

  public :: Mesh3d_Metrics

  !-----------------------------------------------------------------------------
  !> 3D mesh metrics
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
  !>
  !> ### Components referring to points 'l,m' on face `f` of element `e`
  !>
  !>   -  `a      (l,m,f,e)  `  : area coefficient a =|∂x/∂τ₁ × ∂x/∂τ₂|
  !>   -  `n      (l,m,f,e,:)`  : unit normal vector n = (∂x/∂τ₁ × ∂x/∂τ₂)/a
  !>   -  `Ji_n_a (l,m,f,e,:)`  : J⁻¹⋅n

  type Mesh3d_Metrics
    class(Mesh3d_Partition), pointer :: mesh     !< related mesh partition
    type(StandardOperators_1D) :: standard_op     !< standard element operators
    real(RNP), allocatable :: x    (:,:,:,:,:)   !< mesh points
    real(RNP), allocatable :: Jm   (:,:,:,:,:,:) !< Jacobian matrix
    real(RNP), allocatable :: Ji   (:,:,:,:,:,:) !< Jacobian matrix inverse
    real(RNP), allocatable :: Jd   (:,:,:,:)     !< Jacobian matrix determinant
    real(RNP), allocatable :: G    (:,:,:,:,:)   !< Laplacian metrics
    real(RNP), allocatable :: a    (:,:,:,:)     !< area coefficient
    real(RNP), allocatable :: n    (:,:,:,:,:)   !< unit normal vector
    real(RNP), allocatable :: Ji_n (:,:,:,:,:)   !< J⁻¹⋅n
  end type Mesh3d_Metrics

  ! constructor
  interface Mesh3d_Metrics
    module procedure New_Mesh3d_Metrics
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Mesh3d_Metrics constructor

  function New_Mesh3d_Metrics(mesh, po, basis) result(this)
    class(Mesh3d_Partition), target, intent(in) :: mesh  !< mesh partition
    integer,                         intent(in) :: po    !< polynomial order
    character(len=*),      optional, intent(in) :: basis !< 'GL' or 'GLL' ['GLL']

    type(Mesh3d_Metrics) :: this

    call Init_Mesh3d_Metrics(this, mesh, po, basis)

  end function New_Mesh3d_Metrics

  !-----------------------------------------------------------------------------
  !> Initialization of mesh metrics

  subroutine Init_Mesh3d_Metrics(this, mesh, po, basis)

    class(Mesh3d_Metrics), intent(inout) :: this

    class(Mesh3d_Partition), target, intent(in) :: mesh  !< mesh partition
    integer,                         intent(in) :: po    !< polynomial order
    character(len=*),      optional, intent(in) :: basis !< 'GL' or 'GLL' ['GLL']

    ! local variables ..........................................................

    real(RNP), allocatable :: Ds_t(:,:), grad_x(:,:,:,:,:)
    real(RNP) :: a, c, Jd, Jm(3,3), Ji(3,3), n(3)
    integer   :: e, f, i, j, k, p, q, np

    ! basic initialization .....................................................

    this % mesh => mesh
    this % standard_op = StandardOperators_1D(po, basis)
    call mesh % GetPoints(po, basis, this % x)

    allocate( this % Jm (0:po,0:po,0:po,mesh%n_elem,3,3) )
    allocate( this % Ji (0:po,0:po,0:po,mesh%n_elem,3,3) )
    allocate( this % Jd (0:po,0:po,0:po,mesh%n_elem)     )
    allocate( this % G  (0:po,0:po,0:po,mesh%n_elem,6)   )

    allocate( this % a    (0:po,0:po,6,mesh%n_elem)   )
    allocate( this % n    (0:po,0:po,6,mesh%n_elem,3) )
    allocate( this % Ji_n (0:po,0:po,6,mesh%n_elem,3) )

    ! auxiliary data
    allocate(Ds_t(0:po,0:po), grad_x(0:po,0:po,0:po,3,3))
    np = po + 1

    associate(Ds => this % standard_op % D)

      Ds_t = transpose(Ds)

      do e = 1, mesh % n_elem

        ! grad_x = dx/dξ .......................................................

        call TPO_Grad_S(np, Ds, Ds_t, this % x(:,:,:,e,1), grad_x(:,:,:,:,1))
        call TPO_Grad_S(np, Ds, Ds_t, this % x(:,:,:,e,2), grad_x(:,:,:,:,2))
        call TPO_Grad_S(np, Ds, Ds_t, this % x(:,:,:,e,3), grad_x(:,:,:,:,3))

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

    end associate

    ! finalization .............................................................

    deallocate(grad_x)

  end subroutine Init_Mesh3d_Metrics

  !-----------------------------------------------------------------------------
  !> Standard element gradient -- to be replaced by optimized TPO routine

  subroutine TPO_Grad_S(np, Ds, Ds_t, u, v)
    integer,   intent(in)  :: np             !< num points per direction
    real(RNP), intent(in)  :: Ds  (np, np)   !< standard diff matrix
    real(RNP), intent(in)  :: Ds_t(np, np)   !< transposed standard diff matrix
    real(RNP), intent(in)  :: u (np,np,np)   !< scalar element variable
    real(RNP), intent(out) :: v (np,np,np,3) !< gradient of u

    integer :: i, j, k, p

    ! v1 = du/dx1, v2 = du/dx2 .................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,1) = 0
      v(i,j,k,2) = 0
      do p = 1, np
        v(i,j,k,1) = v(i,j,k,1) + Ds_t(p,i) * u(p,j,k)
        v(i,j,k,2) = v(i,j,k,2) + Ds_t(p,j) * u(i,p,k)
      end do
    end do
    end do
    end do

    ! v3 = du/dx3 ..............................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,3) = 0
      do p = 1, np
        v(i,j,k,3) =  v(i,j,k,3) + Ds_t(p,j) * u(i,j,p)
      end do
    end do
    end do
    end do

  end subroutine TPO_Grad_S

  !=============================================================================

end module Mesh_3d__Metrics
