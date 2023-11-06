!> summary:  Application of velocity boundary conditions
!> author:   Joerg Stiller
!> date:     2022/12/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyVelocityBC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of velocity boundary conditions to trace variables
  !>
  !> Construction of outer velocity traces and viscous flux boundary conditions.
  !> The velocity traces `v⁺` are computed using the inner traces `tr_v = v⁻`
  !> given on input, and the boundary values of `v` defined in `bv_u`. In the
  !> absence of the latter, homogeneous conditions are assumed.
  !>
  !> The trace variable `tr_s` is used to return boundary conditions for the
  !> stress vector at outflow faces. These conditions are computed using the
  !> flow variables given in `bv_u`. If the latter absent, the stress vector
  !> is set to zero.

  module subroutine ApplyVelocityBC(this, tr_v, tr_s, bv_u, extrapolate)
    class(INS_Operator_3D), intent(in) :: this
    !< Navier-Stokes operator
    real(RNP), contiguous, intent(inout) :: tr_v(:,:,:,:,:)
    !< velocity traces `tr_v(np,np,6,ne,3)`
    !!   - input:  `tr_v = v⁻`  on all faces
    !!   - output: `tr_v = v⁺`  on boundary faces
    real(RNP), contiguous, optional, intent(inout) :: tr_s(:,:,:,:,:)
    !< stress vector on element faces `tr_s(np,np,6,nl,3)`
    !!   - input:  `tr_s = 0  `  on all faces,
    !!   - output: `tr_s = s_b`  on outflow faces
    class(BoundaryVariable_3D), optional, intent(in) :: bv_u(:)
    !< boundary values depending on BC type specified in `this % bc_v`:
    !!   - `'D'`:  `v  ` in components 1:3
    !!   - `'O'`:  `v,p` in components 1:4

    character(len=*), optional, intent(in) :: extrapolate
    !< switch to replace boundary conditions by extrapolation:
    !!   - `D` :  Dirichlet conditions
    !!   - `O` :  outflow conditions
    !!   - `A` :  all

    real(RNP), allocatable :: theta(:,:), vn(:,:), vv(:,:)
    real(RNP) :: v_max = 0 ! shared with OpenMP !
    real(RNP) :: cv
    integer   :: b, c, e, f, m
    logical   :: ex_db, ex_ob

    if (present(extrapolate)) then
      ex_db = scan(extrapolate,'AD') > 0
      ex_ob = scan(extrapolate,'AO') > 0
    else
      ex_db = .false.
      ex_ob = .false.
    end if

    if (present(tr_s)) then
      allocate(theta(size(tr_s,1), size(tr_s,1)))
      allocate(vn, mold = theta)
      allocate(vv, mold = theta)
    end if

    do b = 1, this % mesh % n_bound
      associate(boundary => this % mesh % boundary(b))

        select case(this % bc_v(b))

        case('D')

          ! Dirichlet conditions: tr_v = v⁺ = 2vᵇ - v⁻ .........................

          if (.not. ex_db) then
            !$omp do
            do f = 1, boundary % n_face
              e = boundary % face(f) % element_id
              m = boundary % face(f) % element_face
              if (present(bv_u)) then
                do c = 1, 3
                  tr_v(:,:,m,e,c) = 2 * bv_u(b) % val(:,:,f,c) - tr_v(:,:,m,e,c)
                end do
              else
                do c = 1, 3
                  tr_v(:,:,m,e,c) = -tr_v(:,:,m,e,c)
                end do
              end if
            end do
          end if

        case('O')

          ! outflow conditions: sᵇ = pn + (nv⋅v + n⋅vv)/2 Θ(n⋅v)

          if (present(tr_s) .and. present(bv_u) .and. .not. ex_ob) then

            !$omp master
            v_max = 0
            !$omp end master

            !$omp do reduction(max:v_max)
            do f = 1, boundary % n_face
              e = boundary % face(f) % element_id
              m = boundary % face(f) % element_face
              associate( n1 => this % sem_v % metrics % n(:,:,m,e,1) &
                       , n2 => this % sem_v % metrics % n(:,:,m,e,2) &
                       , n3 => this % sem_v % metrics % n(:,:,m,e,3) &
                       , v1 => bv_u(b) % val(:,:,f,1)                &
                       , v2 => bv_u(b) % val(:,:,f,2)                &
                       , v3 => bv_u(b) % val(:,:,f,3)                )

                v_max = max(v_max, maxval(abs(n1 * v1 + n2 * v2 + n3 * v3)))

              end associate
            end do

            cv = 1 / max(this % delta_outflow * v_max, epsilon(ONE))

            !$omp do
            do f = 1, boundary % n_face
              e = boundary % face(f) % element_id
              m = boundary % face(f) % element_face
              associate( n1 => this % sem_v % metrics % n(:,:,m,e,1) &
                       , n2 => this % sem_v % metrics % n(:,:,m,e,2) &
                       , n3 => this % sem_v % metrics % n(:,:,m,e,3) &
                       , v1 => bv_u(b) % val(:,:,f,1)                &
                       , v2 => bv_u(b) % val(:,:,f,2)                &
                       , v3 => bv_u(b) % val(:,:,f,3)                &
                       , p  => bv_u(b) % val(:,:,f,4)                &
                       , s1 => tr_s(:,:,m,e,1)                       &
                       , s2 => tr_s(:,:,m,e,2)                       &
                       , s3 => tr_s(:,:,m,e,3)                       )

                  vn = n1 * v1 + n2 * v2 + n3 * v3
                  vv = v1 * v1 + v2 * v2 + v3 * v3

                  theta = HALF * (ONE - tanh(cv * vn))

                  s1 = p * n1 + HALF * (vv * n1 + vn * v1) * theta
                  s2 = p * n2 + HALF * (vv * n2 + vn * v2) * theta
                  s3 = p * n3 + HALF * (vv * n3 + vn * v3) * theta

                end associate
              end do

            end if

        end select

      end associate
    end do

  end subroutine ApplyVelocityBC

  !=============================================================================

end submodule MP_ApplyVelocityBC
