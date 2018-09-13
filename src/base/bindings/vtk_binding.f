!> summary:  Provides an interface to the C binding of the VTK XML writer
!> author:   Joerg Stiller
!> date:     2014/05/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Provides an interface to the C binding of the VTK XML writer.
!===============================================================================

module VTK_Binding
  use C_Binding
  implicit none
  private

  !=============================================================================
  ! public parameters

  !-----------------------------------------------------------------------------
  ! Kind parameters

  !> integer kind matching vtkIdType
  integer, parameter, public :: C_VTK_ID = C_LONG

  !-----------------------------------------------------------------------------
  ! VTK data objects types

  !> VTK data object type for unstructured grids
  integer(C_INT), parameter, public :: VTK_UNSTRUCTURED_GRID = 4

  !-----------------------------------------------------------------------------
  ! VTK data mode types

  integer(C_INT), parameter, public :: VTK_ASCII    = 0 !< ASCII data mode
  integer(C_INT), parameter, public :: VTK_BINARY   = 1 !< binary data mode
  integer(C_INT), parameter, public :: VTK_APPENDED = 2 !< appended data mode

  !-----------------------------------------------------------------------------
  ! VTK cell types

  ! linear cell types
  integer(C_INT), parameter, public :: &
    VTK_EMPTY_CELL       =  0, & !< empty cell
    VTK_VERTEX           =  1, & !< vertex
    VTK_POLY_VERTEX      =  2, & !< set of vertices
    VTK_LINE             =  3, & !< straight line
    VTK_POLY_LINE        =  4, & !< piecewise linear path
    VTK_TRIANGLE         =  5, & !< linear triangle
    VTK_TRIANGLE_STRIP   =  6, & !< strip of linear triangles
    VTK_POLYGON          =  7, & !< linear polygon
    VTK_PIXEL            =  8, & !< pixel, square
    VTK_QUAD             =  9, & !< linear quadrangle
    VTK_TETRA            = 10, & !< linear tetrahedron
    VTK_VOXEL            = 11, & !< voxel, rectangular hexahedron, brick
    VTK_HEXAHEDRON       = 12, & !< linear hexahedron
    VTK_WEDGE            = 13, & !< linear wedge
    VTK_PYRAMID          = 14, & !< linear pyramid
    VTK_PENTAGONAL_PRISM = 15, & !< linear pentagonal prism
    VTK_HEXAGONAL_PRISM  = 16    !< linear hexagonal prism

  ! VTK quadratic isoparametric cell types
  integer(C_INT), parameter, public :: &
    VTK_QUADRATIC_EDGE                   = 21, & !< quadratic edge
    VTK_QUADRATIC_TRIANGLE               = 22, & !< quadratic triangle
    VTK_QUADRATIC_QUAD                   = 23, & !< quadratic quadrangle
    VTK_QUADRATIC_TETRA                  = 24, & !< quadratic tetra
    VTK_QUADRATIC_HEXAHEDRON             = 25, & !< quadratic hexahedron
    VTK_QUADRATIC_WEDGE                  = 26, & !< quadratic wedge
    VTK_QUADRATIC_PYRAMID                = 27, & !< quadratic pyramid
    VTK_BIQUADRATIC_QUAD                 = 28, & !< biquadratic quadrangle
    VTK_TRIQUADRATIC_HEXAHEDRON          = 29, & !< triquadratic hexahedron
    VTK_QUADRATIC_LINEAR_QUAD            = 30, & !< quadratic linear quadrangle
    VTK_QUADRATIC_LINEAR_WEDGE           = 31, & !< quadratic linear wedge
    VTK_BIQUADRATIC_QUADRATIC_WEDGE      = 32, & !< biquadratic quadratic wedge
    VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33, & !< biquadratic quad. hexahedron
    VTK_BIQUADRATIC_TRIANGLE             = 34    !< biquadratic triangle

  !=============================================================================
  ! private parameters

  !-----------------------------------------------------------------------------
  ! VTK data types

  integer(C_INT), parameter :: VTK_INT    =  6
  integer(C_INT), parameter :: VTK_FLOAT  = 10
  integer(C_INT), parameter :: VTK_DOUBLE = 11

  !=============================================================================
  ! public procedures

  public :: VTK_XMLWriter_New
  public :: VTK_XMLWriter_Delete
  public :: VTK_XMLWriter_SetDataObjectType
  public :: VTK_XMLWriter_SetDataModeType
  public :: VTK_XMLWriter_SetExtent
  public :: VTK_XMLWriter_SetPoints
  public :: VTK_XMLWriter_SetOrigin
  public :: VTK_XMLWriter_SetSpacing
  public :: VTK_XMLWriter_SetCoordinates
  public :: VTK_XMLWriter_SetCellsWithType
  public :: VTK_XMLWriter_SetCellsWithTypes
  public :: VTK_XMLWriter_SetPointData
  public :: VTK_XMLWriter_SetCellData
  public :: VTK_XMLWriter_SetFileName
  public :: VTK_XMLWriter_Write
  public :: VTK_XMLWriter_SetNumberOfTimeSteps
  public :: VTK_XMLWriter_Start
  public :: VTK_XMLWriter_WriteNextTimeStep
  public :: VTK_XMLWriter_Stop

  !---------------------------------------------------------------------------
  !> \brief   Create a new instance of a C vtkXMLWriter object.
  !> \author  Joerg Stiller

  interface
    subroutine VTK_XMLWriter_New(writer) &
        bind(C, name='vtkXMLWriterF_New')
      import :: C_INT
      integer(C_INT), intent(inout) :: writer !< writer object
    end subroutine VTK_XMLWriter_New
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Delete the writer object.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This should not be called between VTK_XMLWriter_Start and
  !> VTK_XMLWriter_Stop calls.

  interface
    subroutine VTK_XMLWriter_Delete(writer) &
        bind(C, name='vtkXMLWriterF_Delete')
      import :: C_INT
      integer(C_INT), intent(inout) :: writer !< writer object
    end subroutine VTK_XMLWriter_Delete
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Set the VTK data object type that will be written.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This initializes an empty data object of the given type.
  !> This must be set before setting geometry or data information
  !> can be set only once per writer object.

  interface
    subroutine VTK_XMLWriter_SetDataObjectType(writer, objType) &
        bind(C, name='vtkXMLWriterF_SetDataObjectType')
      import :: C_INT
      integer(C_INT), intent(in) :: writer  !< writer object
      integer(C_INT), intent(in) :: objType !< VTK data object type
    end subroutine VTK_XMLWriter_SetDataObjectType
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Set the VTK writer data mode.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> The following data modes can be specified:
  !>   *  VTK_ASCII
  !>   *  VTK_BINARY
  !>   *  VTK_APPENDED, which is the default

  interface
    subroutine VTK_XMLWriter_SetDataModeType(writer, dataModeType) &
        bind(C, name='vtkXMLWriterF_SetDataModeType')
      import :: C_INT
      integer(C_INT), intent(in) :: writer       !< writer object
      integer(C_INT), intent(in) :: dataModeType !< data mode
    end subroutine VTK_XMLWriter_SetDataModeType
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Set the extent of a structured data set.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may be used only after SetDataObjectType has been called
  !> with a structured data object type.

  interface
    subroutine VTK_XMLWriter_SetExtent(writer, extent) &
        bind(C, name='vtkXMLWriterF_SetExtent')
      import :: C_INT
      integer(C_INT), intent(in) :: writer    !< writer object
      integer(C_INT), intent(in) :: extent(6) !< data set extent
    end subroutine VTK_XMLWriter_SetExtent
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Set the points of a point data set.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> May not be used for data objects with implicit points.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface VTK_XMLWriter_SetPoints
    module procedure VTK_XMLWriter_SetPoints_F ! C_FLOAT
    module procedure VTK_XMLWriter_SetPoints_D ! C_DOUBLE
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Set the origin of an image data set.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be used for image data.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface
    subroutine VTK_XMLWriter_SetOrigin(writer, origin) &
        bind(C, name='vtkXMLWriterF_SetOrigin')
      import :: C_INT, C_DOUBLE
      integer(C_INT), intent(in) :: writer    !< writer object
      real(C_DOUBLE), intent(in) :: origin(3) !< Cart. coordinates of origin
    end subroutine VTK_XMLWriter_SetOrigin
  end interface

  !---------------------------------------------------------------------------
  !> \brief   Set the spacing of an image data set.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be used for image data.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface
    subroutine VTK_XMLWriter_SetSpacing(writer, spacing) &
        bind(C, name='vtkXMLWriterF_SetSpacing')
      import :: C_INT, C_DOUBLE
      integer(C_INT), intent(in) :: writer     !< writer object
      real(C_DOUBLE), intent(in) :: spacing(3) !< spacing of the image data set
    end subroutine VTK_XMLWriter_SetSpacing
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Set the coordinates along one axis of a rectilinear grid.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be used for rectilinear grids.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface VTK_XMLWriter_SetCoordinates
    module procedure VTK_XMLWriter_SetCoordinates_F ! C_FLOAT
    module procedure VTK_XMLWriter_SetCoordinates_D ! C_DOUBLE
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Set a point data array by name.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> The name of the array is required and should describe the purpose
  !> of the data. Use dataType to specify the scalar type used in the
  !> given data array. Use numTuples to specify the number of tuples
  !> and numComponents to specify the number of scalar components in
  !> each tuple.
  !>
  !> The data array must be shaped [numComponents,numTuples] where numTuples
  !> equals the number of points indicated by SetExtent and/or SetPoints.
  !>
  !> The role can be one of "SCALARS", "VECTORS", "NORMALS", "TENSORS",
  !> or "TCOORDS" and specifies that the array should be designated as
  !> the active array for the named role. Other values for role are
  !> ignored.
  !>
  !> This may be used for all data types.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface VTK_XMLWriter_SetPointData
    module procedure VTK_XMLWriter_SetPointData_I1 ! C_INT, one component
    module procedure VTK_XMLWriter_SetPointData_I2 ! C_INT, multicomponent
    module procedure VTK_XMLWriter_SetPointData_F1 ! C_FLOAT, one component
    module procedure VTK_XMLWriter_SetPointData_F2 ! C_FLOAT, multicomponent
    module procedure VTK_XMLWriter_SetPointData_D1 ! C_DOUBLE, one component
    module procedure VTK_XMLWriter_SetPointData_D2 ! C_DOUBLE, multicomponent
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Set a cell data array by name.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> The name of the array is required and should describe the purpose
  !> of the data. Use dataType to specify the scalar type used in the
  !> given data array. Use numTuples to specify the number of tuples
  !> and numComponents to specify the number of scalar components in
  !> each tuple.
  !>
  !> The data array must be shaped [numComponents,numTuples] where numTuples
  !> equals the number of cells set by SetCells.
  !>
  !> The role can be one of "SCALARS", "VECTORS", "NORMALS", "TENSORS",
  !> or "TCOORDS" and specifies that the array should be designated as
  !> the active array for the named role. Other values for role are
  !> ignored.
  !>
  !> This may be used for all data types.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface VTK_XMLWriter_SetCellData
    module procedure VTK_XMLWriter_SetCellData_I1 ! C_INT, one component
    module procedure VTK_XMLWriter_SetCellData_I2 ! C_INT, multicomponent
    module procedure VTK_XMLWriter_SetCellData_F1 ! C_FLOAT, one component
    module procedure VTK_XMLWriter_SetCellData_F2 ! C_FLOAT, multicomponent
    module procedure VTK_XMLWriter_SetCellData_D1 ! C_DOUBLE, one component
    module procedure VTK_XMLWriter_SetCellData_D2 ! C_DOUBLE, multicomponent
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Write the data to a file immediately. This is not used when
  !>          writing time-series data. Returns 1 for success and 0 for failure.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be called after SetFileName and SetDataObjectType.

  interface
    subroutine VTK_XMLWriter_Write(writer, success) &
        bind(C, name='vtkXMLWriterF_Write')
      import :: C_INT
      integer(C_INT), intent(in)  :: writer
      integer(C_INT), intent(out) :: success
    end subroutine VTK_XMLWriter_Write
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Sets the number of time steps that will be written.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> Sets the number of time steps that will be written between upcoming
  !> VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.
  !> This is used when writing time-series data.
  !> This may be used for all data types.
  !> It may not be called before VTK_XMLWriter_SetDataObjectType
  !> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

  interface
    subroutine VTK_XMLWriter_SetNumberOfTimeSteps(writer, numTimeSteps) &
        bind(C, name='vtkXMLWriterF_SetNumberOfTimeSteps')
      import :: C_INT
      integer(C_INT), intent(in) :: writer
      integer(C_INT), intent(in) :: numTimeSteps
    end subroutine VTK_XMLWriter_SetNumberOfTimeSteps
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Start writing a time-series to the output file.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be called after SetFileName, SetDataObjectType, and
  !> SetNumberOfTimeSteps. It may not be called a second time until
  !> after an intervening call to VTK_XMLWriter_Stop.

  interface
    subroutine VTK_XMLWriter_Start(writer) &
        bind(C, name='vtkXMLWriterF_Start')
      import :: C_INT
      integer(C_INT), intent(in) :: writer
    end subroutine VTK_XMLWriter_Start
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Write one time step of a time-series to the output file.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> The current data set by SetPointData and SetCellData will be written.
  !> Use timeValue to specify the time associated with the time step
  !> being written. \n
  !> This may only be called after VTK_XMLWriter_Start has been called.
  !> It should be called NumberOfTimeSteps times before calling
  !> VTK_XMLWriter_Stop.

  interface
    subroutine VTK_XMLWriter_WriteNextTimeStep(writer, timeValue) &
        bind(C, name='vtkXMLWriterF_WriteNextTimeStep')
      import :: C_INT, C_DOUBLE
      integer(C_INT), intent(in) :: writer
      real(C_DOUBLE), intent(in) :: timeValue
    end subroutine VTK_XMLWriter_WriteNextTimeStep
  end interface

  !-----------------------------------------------------------------------------
  !> \brief   Stop writing a time-series to the output file.
  !> \author  Joerg Stiller
  !>
  !> \details
  !> This may only be called after VTK_XMLWriter_Start and NumberOfTimeSteps
  !> calls to VTK_XMLWriter_WriteNextTimeStep.

  interface
    subroutine VTK_XMLWriter_Stop(writer) &
        bind(C, name='vtkXMLWriterF_Stop')
      import :: C_INT
      integer(C_INT), intent(in) :: writer
    end subroutine VTK_XMLWriter_Stop
  end interface

contains

!===============================================================================
! Specific implementations of VTK_XMLWriter_SetPoints

!-------------------------------------------------------------------------------
!> \brief   Set the points of a real point data set with kind C_FLOAT.
!> \author  Joerg Stiller
!>
!> \details
!> The array x(1:3,:) specifies the Cartesian coordinates of the 3D points.

subroutine VTK_XMLWriter_SetPoints_F(writer, x)
  integer(C_INT), intent(in) :: writer !< writer object
  real(C_FLOAT),  intent(in) :: x(:,:) !< Cartesian point coordinates

  interface
    subroutine SetPoints(writer, dataType, x, numPoints) &
        bind(C, name='vtkXMLWriterF_SetPoints')
      import :: C_INT, C_VTK_ID, C_FLOAT
      integer(C_INT),    intent(in) :: writer     !< writer object
      integer(C_INT),    intent(in) :: dataType   !< VTK data type
      real(C_FLOAT),     intent(in) :: x(3,*)     !< coordinates
      integer(C_VTK_ID), intent(in) :: numPoints  !< number of points
    end subroutine SetPoints
  end interface

  call SetPoints(writer, VTK_FLOAT, x, size(x,2,C_VTK_ID))

end subroutine VTK_XMLWriter_SetPoints_F

!-------------------------------------------------------------------------------
!> \brief   Set the points of a real point data set with kind C_DOUBLE
!> \author  Joerg Stiller
!>
!> \details
!> The array x(1:3,:) specifies the Cartesian coordinates of the 3D points.

subroutine VTK_XMLWriter_SetPoints_D(writer, x)
  integer(C_INT), intent(in) :: writer !< writer object
  real(C_DOUBLE), intent(in) :: x(:,:) !< Cartesian point coordinates

  interface
    subroutine SetPoints(writer, dataType, x, numPoints) &
        bind(C, name='vtkXMLWriterF_SetPoints')
      import :: C_INT, C_VTK_ID, C_DOUBLE
      integer(C_INT),    intent(in) :: writer    !< writer object
      integer(C_INT),    intent(in) :: dataType  !< VTK data type
      real(C_DOUBLE),    intent(in) :: x(3,*)    !< coordinates
      integer(C_VTK_ID), intent(in) :: numPoints !< number of points
    end subroutine SetPoints
  end interface

  call SetPoints(writer, VTK_DOUBLE, x, size(x,2,C_VTK_ID))

end subroutine VTK_XMLWriter_SetPoints_D

!=============================================================================
! Specific implementations of VTK_XMLWriter_SetCoordinates

!-------------------------------------------------------------------------------
!> \brief   Set the coordinates along one axis of a rectilinear grid data set.
!> \author  Joerg Stiller
!>
!> \details
!> Specify axis 0 for X, 1 for Y, and 2 for Z. Use dataType to specify
!> the scalar type used in the given data array.

subroutine VTK_XMLWriter_SetCoordinates_F(writer, axis, data)
  integer(C_INT), intent(in) :: writer  !< writer object
  integer(C_INT), intent(in) :: axis    !< axis identifier
  real(C_FLOAT),  intent(in) :: data(:) !< grid points along the axis

  interface
    subroutine SetCoordinates(writer, axis, dataType, data, numCoordinates) &
        bind(C, name='vtkXMLWriterF_SetCoordinates')
      import :: C_INT, C_VTK_ID, C_FLOAT
      integer(C_INT),    intent(in) :: writer
      integer(C_INT),    intent(in) :: axis
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numCoordinates
      real(C_FLOAT),     intent(in) :: data(numCoordinates)
    end subroutine SetCoordinates
  end interface

  call SetCoordinates(writer, axis, VTK_FLOAT, data, size(data,1,C_VTK_ID))

end subroutine VTK_XMLWriter_SetCoordinates_F

!-------------------------------------------------------------------------------
!> \brief   Set the coordinates along one axis of a rectilinear grid data set.
!> \author  Joerg Stiller
!>
!> \details
!> Specify axis 0 for X, 1 for Y, and 2 for Z. Use dataType to specify the
!> scalar type used in the given data array.

subroutine VTK_XMLWriter_SetCoordinates_D(writer, axis, data)
  integer(C_INT), intent(in) :: writer  !< writer object
  integer(C_INT), intent(in) :: axis    !< axis identifier
  real(C_DOUBLE), intent(in) :: data(:) !< grid points along the axis

  interface
    subroutine SetCoordinates(writer, axis, dataType, data, &
        numCoordinates) bind(C, name='vtkXMLWriterF_SetCoordinates')
      import :: C_INT, C_VTK_ID, C_DOUBLE
      integer(C_INT),    intent(in) :: writer
      integer(C_INT),    intent(in) :: axis
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numCoordinates
      real(C_DOUBLE),    intent(in) :: data(numCoordinates)
    end subroutine SetCoordinates
  end interface

  call SetCoordinates(writer, axis, VTK_DOUBLE, data, size(data,1,C_VTK_ID))

end subroutine VTK_XMLWriter_SetCoordinates_D

!=============================================================================

!-----------------------------------------------------------------------------
!> \brief   Set a cell array on the data object. All cells have the same type.
!> \author  Joerg Stiller
!>
!> \details
!> For unstructured grid data objects, the cellType can be any type.
!> For polygonal data objects, the cellType must be
!>   * VTK_VERTEX,
!>   * VTK_POLY_VERTEX,
!>   * VTK_LINE,
!>   * VTK_POLY_LINE,
!>   * VTK_TRIANGLE,
!>   * VTK_TRIANGLE_STRIP, or
!>   * a cyclically connected simple cell type such as VTK_POLYGON.
!>
!> The cell array is dimensioned as [0:M,N], where M is the number of points
!> in the cell and N the number of cells. Each line uses the M+1 entries
!> [M,id1,id2,...,idM] where id1, ..., idM are the point IDs.
!>
!> This may only be used for unstructured grid and polygonal data types.
!> It may not be called before VTK_XMLWriter_SetDataObjectType
!> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

subroutine VTK_XMLWriter_SetCellsWithType(writer, cell_type, cell)
  integer(C_INT),    intent(in) :: writer    !< writer object
  integer(C_INT),    intent(in) :: cell_type !< VTK cell type
  integer(C_VTK_ID), intent(in) :: cell(:,:) !< cell array

  interface
    subroutine SetCellsWithType(writer, cell_type, n_cell, cell, cell_size) &
        bind(C, name='vtkXMLWriterF_SetCellsWithType')
      import :: C_INT, C_VTK_ID
      integer(C_INT),    intent(in) :: writer            !< writer object
      integer(C_INT),    intent(in) :: cell_type         !< cell type
      integer(C_VTK_ID), intent(in) :: cell_size         !< cell size
      integer(C_VTK_ID), intent(in) :: n_cell            !< number of cells
      integer(C_VTK_ID), intent(in) :: cell(cell_size,*) !< cell array
    end subroutine SetCellsWithType
  end interface

  call SetCellsWithType(writer, cell_type, size(cell,2,C_VTK_ID), cell, &
      size(cell,1,C_VTK_ID))

end subroutine VTK_XMLWriter_SetCellsWithType

!-------------------------------------------------------------------------------
!> \brief   Set a cell array on the data object. Each cell has its own type.
!> \author  Joerg Stiller
!>
!> \details
!> The cell array is dimensioned as [0:M,N], where M is the number of points
!> in the cell and N the number of cells. Each line uses the M+1 entries
!> [M,id1,id2,...,idM] where id1, ..., idM are the point IDs.
!>
!> This may only be used for unstructured grid data objects.
!> It may not be called before VTK_XMLWriter_SetDataObjectType
!> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

subroutine VTK_XMLWriter_SetCellsWithTypes(writer, cell_type, cell)
  integer(C_INT),    intent(in) :: writer       !< writer object
  integer(C_INT),    intent(in) :: cell_type(:) !< VTK cell types
  integer(C_VTK_ID), intent(in) :: cell(:,:)    !< cell array

  interface
    subroutine SetCellsWithTypes(writer, cell_type, n_cell, cell, cell_size) &
        bind(C, name='vtkXMLWriterF_SetCellsWithTypes')
      import :: C_INT, C_VTK_ID
      integer(C_INT),    intent(in) :: writer            !< writer object
      integer(C_INT),    intent(in) :: cell_type(*)      !< cell type
      integer(C_VTK_ID), intent(in) :: cell_size         !< cell size
      integer(C_VTK_ID), intent(in) :: n_cell            !< number of cells
      integer(C_VTK_ID), intent(in) :: cell(cell_size,*) !< cell array
    end subroutine SetCellsWithTypes
  end interface

  call SetCellsWithTypes(writer, cell_type, size(cell,2,C_VTK_ID), cell, &
      size(cell,1,C_VTK_ID))

end subroutine VTK_XMLWriter_SetCellsWithTypes

!=============================================================================
! Specific implementations of VTK_XMLWriter_SetPointData

!-------------------------------------------------------------------------------
!> \brief   Set a scalar point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_I1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  integer(C_INT),   intent(in) :: data(:) !< point data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      integer(C_INT),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_INT, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetPointData_I1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_I2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  integer(C_INT),   intent(in) :: data(:,:) !< point data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      integer(C_INT),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_INT, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetPointData_I2

!-------------------------------------------------------------------------------
!> \brief   Set a scalar point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_F1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  real(C_FLOAT),    intent(in) :: data(:) !< point data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_FLOAT
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_FLOAT),     intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_FLOAT, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetPointData_F1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_F2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  real(C_FLOAT),    intent(in) :: data(:,:) !< point data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_FLOAT
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_FLOAT),     intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_FLOAT, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetPointData_F2

!-------------------------------------------------------------------------------
!> \brief   Set a scalar point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_D1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  real(C_DOUBLE),   intent(in) :: data(:) !< point data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_DOUBLE
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_DOUBLE),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_DOUBLE, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetPointData_D1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component point data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetPointData_D2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  real(C_DOUBLE),   intent(in) :: data(:,:) !< point data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetPointData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetPointData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_DOUBLE
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_DOUBLE),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetPointData
  end interface

  call SetPointData(writer, CString(name), VTK_DOUBLE, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetPointData_D2

!=============================================================================
! Specific implementations of VTK_XMLWriter_SetCellData

!-------------------------------------------------------------------------------
!> \brief   Set a scalar cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_I1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  integer(C_INT),   intent(in) :: data(:) !< cell data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      integer(C_INT),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_INT, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetCellData_I1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_I2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  integer(C_INT),   intent(in) :: data(:,:) !< cell data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      integer(C_INT),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_INT, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetCellData_I2

!-------------------------------------------------------------------------------
!> \brief   Set a scalar cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_F1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  real(C_FLOAT),    intent(in) :: data(:) !< cell data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_FLOAT
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_FLOAT),     intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_FLOAT, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetCellData_F1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_F2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  real(C_FLOAT),    intent(in) :: data(:,:) !< cell data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_FLOAT
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_FLOAT),     intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_FLOAT, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetCellData_F2

!-------------------------------------------------------------------------------
!> \brief   Set a scalar cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_D1(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer  !< writer object
  character(len=*), intent(in) :: name    !< name of the data set
  real(C_DOUBLE),   intent(in) :: data(:) !< cell data set
  character(len=*), intent(in) :: role    !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_DOUBLE
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_DOUBLE),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_DOUBLE, data, &
    size(data,1,C_VTK_ID), 1_C_INT, CString(role))

end subroutine VTK_XMLWriter_SetCellData_D1

!-------------------------------------------------------------------------------
!> \brief   Set a multi-component cell data array by name.
!> \author  Joerg Stiller

subroutine VTK_XMLWriter_SetCellData_D2(writer, name, data, role)
  integer(C_INT),   intent(in) :: writer    !< writer object
  character(len=*), intent(in) :: name      !< name of the data set
  real(C_DOUBLE),   intent(in) :: data(:,:) !< cell data set
  character(len=*), intent(in) :: role      !< role of the data set

  interface
    subroutine SetCellData(writer, name, dataType, data, numTuples, &
        numComponents, role) bind(C, name='vtkXMLWriterF_SetCellData')
      import :: C_INT, C_VTK_ID, C_CHAR, C_DOUBLE
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT),    intent(in) :: dataType
      integer(C_VTK_ID), intent(in) :: numTuples
      integer(C_INT),    intent(in) :: numComponents
      real(C_DOUBLE),    intent(in) :: data(numComponents,numTuples)
      character(C_CHAR), intent(in) :: role(*)
    end subroutine SetCellData
  end interface

  call SetCellData(writer, CString(name), VTK_DOUBLE, data, &
    size(data,2,C_VTK_ID), size(data,1,C_INT), CString(role))

end subroutine VTK_XMLWriter_SetCellData_D2

!=============================================================================

!-----------------------------------------------------------------------------
!> \brief   Set the name of the file into which the data are to be written.
!> \author  Joerg Stiller
!>
!> \details
!> This may be used for all data types.
!> It may not be called before VTK_XMLWriter_SetDataObjectType
!> between VTK_XMLWriter_Start and VTK_XMLWriter_Stop calls.

subroutine VTK_XMLWriter_SetFileName(writer, fileName)
  integer(C_INT),   intent(in) :: writer   !< writer object
  character(len=*), intent(in) :: fileName !< file name

  interface
    subroutine SetFileName(writer, fileName) &
        bind(C, name='vtkXMLWriterF_SetFileName')
      import :: C_INT, C_CHAR
      integer(C_INT),    intent(in) :: writer
      character(C_CHAR), intent(in) :: fileName(*)
    end subroutine SetFileName
  end interface

  call SetFileName(writer, CString(fileName))

end subroutine VTK_XMLWriter_SetFileName

!=============================================================================

end module VTK_Binding