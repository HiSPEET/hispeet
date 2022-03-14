  !=============================================================================
  ! External module procedures

  interface

    !---------------------------------------------------------------------------
    !> Restrict a structured mesh variable `v` to subdomains

    module subroutine Restrict_S(this, mesh, buf_v, v, vs)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_v
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< extended mesh variable
      real(RWP),      intent(out)   :: vs(:,:,:,:)   !< restricted variable
    end subroutine Restrict_S

    !---------------------------------------------------------------------------
    !> Restrict an unstructured mesh variable `v` to subdomains

    module subroutine Restrict_U(this, mesh, buf_v, v, vs)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_v
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< extended mesh variable
      real(RWP),      intent(out)   :: vs(:,:,:,:)   !< restricted variable
    end subroutine Restrict_U

    !-----------------------------------------------------------------------------
    !>  Merge of subdomain core contributions into mesh variable

    module subroutine Merge_C(this, mesh, vs, v)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< subdomain solutions
      real(RWP),      intent(inout) :: vs(:,:,:,:)   !< mesh variable
    end subroutine Merge_C

    !-----------------------------------------------------------------------------
    !> Merge subdomain contributions into structured mesh variable

    module subroutine Merge_S(this, mesh, buf_vs, vs, v)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_vs
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< subdomain solutions
      real(RWP),      intent(inout) :: vs(:,:,:,:)   !< mesh variable
    end subroutine Merge_S

    !-----------------------------------------------------------------------------
    !> Merge subdomain contributions into unstructured mesh variable

    module subroutine Merge_U(this, mesh, buf_vs, vs, v)
      class(DG_SchwarzOperator_3D), intent(in) :: this
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_vs
      class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
      real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< subdomain solutions
      real(RWP),      intent(inout) :: vs(:,:,:,:)   !< mesh variable
    end subroutine Merge_U

  end interface

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Restrict mesh variable to subdomains
  !>
  !> The mesh variable must be dimensioned as `v(0:po,0:po,0:po,ne+ng)`, where
  !> `ne` is the number of local elements and `ng` the number of ghosts.

  subroutine RestrictResidual_RWP(this, mesh, buf_v, v, vs)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_v
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: vs(:,:,:,:)   !< restricted variable

    if (mesh % structured) then
      call Restrict_S(this, mesh, buf_v, v, vs)
    else
      call Restrict_U(this, mesh, buf_v, v, vs)
    end if

  end subroutine RestrictResidual_RWP

  !-----------------------------------------------------------------------------
  !> Merge subdomain contributions into mesh variable
  !>
  !> The subdomain variable must be dimensioned as `v(ns,ns,ns,ne+ng)`, where
  !> `ns` is the numbeer of subdomain points per direction, `ne` the number of
  !> local elements and `ng` the number of ghosts.

  subroutine MergeCorrections_RWP(this, mesh, buf_vs, vs, v)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_vs
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< subdomain solutions
    real(RWP),      intent(inout) :: vs(:,:,:,:)   !< mesh variable

    if (this % restrictive) then ! merge core regions only
      call Merge_C(this, mesh, vs, v)
    else if (mesh % structured) then
      call Merge_S(this, mesh, buf_vs, vs, v)
    else
      call Merge_U(this, mesh, buf_vs, vs, v)
    end if

  end subroutine MergeCorrections_RWP

  !=============================================================================
