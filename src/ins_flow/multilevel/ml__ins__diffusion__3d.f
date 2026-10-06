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

!> summary:  Multilevel Navier-Stokes diffusion step
!> author:   Joerg Stiller
!> date:     2026/09/30
!===============================================================================

module ML__INS__Diffusion__3D

  use Kind_Parameters
  use Constants
  use XMPI
  use Logging_Levels
  use Array_Reductions

  use Boundary_Variable__3D
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D

  use INS__Operator__3D

  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__Array_Reductions__3D

  use ML__INS__Operator__3D

  implicit none
  private

  public :: ML_INS_Diffusion_3D
  public :: ML_INS_DiffusionOptions_3D

  !-----------------------------------------------------------------------------
  !> Multilevel diffusion for incompressible Navier-Stokes problems

  type ML_INS_Diffusion_3D

    class(ML_INS_Operator_3D), pointer :: ml_ins

    character :: fc_projection   !< fine-to-coarse projection method
    logical   :: fc_reset_twigs  !< reinitialize coarse solution
    integer   :: i_crs           !< max num coarse solver iterations
    integer   :: i_max           !< max number of multigrid iterations
    integer   :: ns_1            !< num pre-smoothing steps
    integer   :: ns_2            !< num post-smoothing steps
    integer   :: ns_c            !< num continuation smoothing steps
    real(RNP) :: r_crs           !< min coarse residual reduction,  if > 0
    real(RNP) :: r_red           !< min residual reduction,  if > 0
    real(RNP) :: r_max           !< max admissible residual, if > 0

  contains

    procedure :: Init_ML_INS_Diffusion_3D
    procedure :: DiffusionResidual
    procedure :: DiffusionStep

  end type ML_INS_Diffusion_3D

  ! constructor
  interface ML_INS_Diffusion_3D
    procedure New_ML_INS_Diffusion_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Multilevel diffusion options

  type ML_INS_DiffusionOptions_3D

    character :: fc_projection  = 'I'     !< projection method {'I','P'}
    logical   :: fc_reset_twigs = .false. !< reset twigs in downward leg

    integer   :: i_crs   = 20 !< max num coarse solver iterations
    integer   :: i_max   =  1 !< max num multigrid iterations (cycles)
    integer   :: ns_1    =  1 !< num pre-smoothing steps
    integer   :: ns_2    =  1 !< num post-smoothing steps
    integer   :: ns_c    = -1 !< num continuation smoothing steps         [auto]
    real(RNP) :: r_crs   = -1 !< min coarse residual reduction, overrides `r_red`
    real(RNP) :: r_red   = -1 !< min residual reduction,  if > 0
    real(RNP) :: r_max   = -1 !< max admissible residual, if > 0

  contains
    procedure :: Bcast => Bcast_ML_INS_DiffusionOptions_3D
  end type ML_INS_DiffusionOptions_3D

  !=============================================================================
  ! Interfaces to submodule procedures

  interface

    !---------------------------------------------------------------------------
    !> FAS residual of the implicit viscous subproblem

    module subroutine DiffusionResidual(this, tau, mu, nu, bv, f, u, r, l_top)
      class(ML_INS_Diffusion_3D),    intent(in)    :: this
      real(RNP),                     intent(in)    :: tau
      class(ML_MeshVariable_3D),     intent(in)    :: mu
      class(ML_MeshVariable_3D),     intent(in)    :: nu
      class(ML_BoundaryVariable_3D), intent(in)    :: bv
      class(ML_MeshVariable_3D),     intent(in)    :: f
      class(ML_MeshVariable_3D),     intent(in)    :: u
      class(ML_MeshVariable_3D),     intent(inout) :: r
      integer,             optional, intent(in)    :: l_top
    end subroutine DiffusionResidual

  !-----------------------------------------------------------------------------
  !> Solution of the implicit viscous subproblem

    module subroutine DiffusionStep(this, tau, mu, nu, bv, f, u, l_top)
      class(ML_INS_Diffusion_3D),    intent(in)    :: this
      real(RNP),                     intent(in)    :: tau
      class(ML_MeshVariable_3D),     intent(in)    :: mu
      class(ML_MeshVariable_3D),     intent(in)    :: nu
      class(ML_BoundaryVariable_3D), intent(in)    :: bv
      class(ML_MeshVariable_3D),     intent(inout) :: f
      class(ML_MeshVariable_3D),     intent(inout) :: u
      integer,             optional, intent(in)    :: l_top
    end subroutine DiffusionStep

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of ML_INS_Diffusion_3D

  !-----------------------------------------------------------------------------
  !> Constructor of ML_INS_Diffusion_3D

  function New_ML_INS_Diffusion_3D(opt, ml_ins) result(this)
    class(ML_INS_DiffusionOptions_3D), intent(in) :: opt
    class(ML_INS_Operator_3D),         intent(in) :: ml_ins
    type(ML_INS_Diffusion_3D) :: this

    call Init_ML_INS_Diffusion_3D(this, opt, ml_ins)

  end function New_ML_INS_Diffusion_3D

  !-----------------------------------------------------------------------------
  !> Initialization of ML_INS_Diffusion_3D

  subroutine Init_ML_INS_Diffusion_3D(this, opt, ml_ins)
    class(ML_INS_Diffusion_3D),        intent(inout) :: this
    class(ML_INS_DiffusionOptions_3D), intent(in)    :: opt
    class(ML_INS_Operator_3D), target, intent(in)    :: ml_ins

    this % ml_ins => ml_ins

    this % fc_projection  = opt % fc_projection
    this % fc_reset_twigs = opt % fc_reset_twigs

    this % i_crs = opt % i_crs
    this % i_max = opt % i_max
    this % ns_1  = opt % ns_1
    this % ns_2  = opt % ns_2

    if (opt % ns_c < 0) then
      this % ns_c = max(this%ns_1, this%ns_2)
    else
      this % ns_c = opt % ns_c
    end if

    this % r_crs = opt % r_crs
    this % r_red = opt % r_red
    this % r_max = opt % r_max

    if (this % r_crs < 0) then
      this % r_crs = this % r_red
    end if

  end subroutine Init_ML_INS_Diffusion_3D

  !=============================================================================
  ! Type-bound procedures of ML_INS_DiffusionOptions_3D

  subroutine Bcast_ML_INS_DiffusionOptions_3D(this, root, comm)
   class(ML_INS_DiffusionOptions_3D), intent(inout) :: this
   integer,        intent(in) :: root !< rank of broadcast root
   type(MPI_Comm), intent(in) :: comm !< MPI communicator

   call XMPI_Bcast(this % fc_projection  , root, comm)
   call XMPI_Bcast(this % fc_reset_twigs , root, comm)
   call XMPI_Bcast(this % i_crs          , root, comm)
   call XMPI_Bcast(this % i_max          , root, comm)
   call XMPI_Bcast(this % ns_1           , root, comm)
   call XMPI_Bcast(this % ns_2           , root, comm)
   call XMPI_Bcast(this % ns_c           , root, comm)
   call XMPI_Bcast(this % r_crs          , root, comm)
   call XMPI_Bcast(this % r_red          , root, comm)
   call XMPI_Bcast(this % r_max          , root, comm)

  end subroutine Bcast_ML_INS_DiffusionOptions_3D

  !=============================================================================

end module ML__INS__Diffusion__3D
