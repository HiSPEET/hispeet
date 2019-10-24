# HiSPEET Style Guide

> _Die Schönheit ist die richtige Übereinstimmung der Teile miteinander und mit dem Ganzen._  
> __Werner Heisenberg__

## Introduction

This guide summarizes the rules for designing for program units in __HiSPEET__, including inline documentation using  [FORD](https://github.com/Fortran-FOSS-Programmers/ford/wiki). As only the principal aspects are covered here, you are strongly encouraged to have a closer look at _real_ modules in `src/base` and example programs in `app/examples`.

## General

### Nomen est omen

__Names__ should transport a meaning and give orientation. Therefore, take your time to find concise but descriptive names for program entities. Additionally, use the  case and underscores to improve readability and to hint at the nature of the entity . 

These are the basic rules:

  - Use lower case Fortran keywords and intrinsic procedures 

  - Program units and types are set in camel case, e.g.
     * `module Gauss_Jacobi`
     * `program FlowSolver`
     * `subroutine SetProblemParameters`
     * `function MeanValue`
 
  - In module names words are separated by underscores
 
  - Subroutine, function and type names are preferably written in camel case with
    no underscores in-between. However, underscores can be inserted for clarity or readability, e.g.
      * `subroutine VTK_Export`
 
  - Constants are written in upper case such as 
      * `PI` for the number $\pi$, or `ONE` for 1, and
      * `RSP`,`RDP`, `RHP`, `RNP` for the real kind parameters defined in module `Kind_Parameters`
      
  - Variables are generally written in lower case with underscores as separators.  However, in favor of conformity to mathematical notation capitalized names can 
be used where appropriate, e.g.
      * `M`, `D`, ` L`  for mass, differentiation and stiffness matrices
      * `i`, `j`, `k`, `l ` etc. for loop indices
      * `c`, `dx`, `v_rms`, `last_exit` 

### Source files

  - Place every program unit in a separate file.

  - Choose the file name to match the contained program unit written in lowercase.

  - Indicate the file type by appending as suffix
      * `.f` for Fortran
      * `.c` and `.h` for C source and header files, respectively
      * `.py` for Python

  - Use the free source format with Fortran.

### Formatting

#### Rules for Fortran sources

  - Improve readability by consistent indentation:
      * The primary program unit, e.g. module or program, always starts in the first column.
      * Statements inside a program unit are indented by two additional columns.
      * An extra indent of two columns is applied to statements inside a block structure except control commands such as `else` or `case`.
      * The extra indent can be skipped for nested `do` loops with no statements placed in-between.
      * Continuation lines are indented with at least two additional columns.
  - Constrain statements and comments to 80 columns.
  - Put only one statement on a line.
  - Use blank spaces to separate variables and operators.
  - Align groups of similar statements.
  - Avoid tabulators.

## Program units

Each program unit is supplemented with a header and inline comments compatible with the [FORD](https://github.com/Fortran-FOSS-Programmers/ford/wiki) documentation generator. _FORD_ comments come in two flavors:

  - The _predoc_ mark `!>` comments the entity defined in the next Fortran.
  - The _docmark_ `!<` starts an inline comment describing the preceding entity
statement.
  - Both comment types can be continued by adding lines starting with the mark or `!!`.

The header of a program, module or submodule is composed of comment lines containing the following fields:

  - `summary:` a short, ideally single-line description
  - `author:`  the author(s)
  - `date:`    creation time in the format yyyy/mm/dd
  - `license:` license information (just copy from existing file)

These fields can be followed be a more detailed description placed in a compound of _predoc_ comment lines.

Optionally special environments such as `@note` or `@todo` can be appended to the header. For more details on the supported Markdown features see the [FORD wiki](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation). 


### Modules

<a name="module_Gauss_Jacobi">_Listing:_</a> Stripped-down version of the `Gauss_Jacobi` module contained in `gauss_jacobi.f`.

```Fortran  
!> summary:  Implementation of Jacobi polynomials and Gauss quadratures
!> author:   Joerg Stiller
!> date:     2013/05/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Implementation of Jacobi polynomials and Gauss quadratures
!>
!> Provides procedures for evaluating
!>
!>   *  the values, derivatives and zeros of Jacobi polynomials,
!>   *  the Gauss-Legendre (GL) points and weights, and the related
!>      Lagrange polynomials and their derivatives,
!>   *  the Gauss-Lobatto-Legendre (GLL) points and weights, and the
!>      related Lagrange polynomials and their derivatives.
!>   *  the left-sided Gauss-Radau-Legendre (GRL) points and weights,
!>      and the related Lagrange polynomials and their derivatives.
!>
!> Implementation follows G.E. Karniadakis & S.J. Sherwin,
!> Spectral/hp Element Methods for CFD. Oxford Univ. Press, 2005.
!==============================================================================

module Gauss_Jacobi
  use Kind_parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE, FOUR, PI
  implicit none
  private

  public :: JacobiPolynomial
  public :: JacobiPolynomialDerivative
  public :: JacobiPolynomialZeros

  ...

  !> Mininmal admissible distance distance to collocation points.
  real(RNP), parameter :: TOL = 1.0e-10_RNP

contains

  !=============================================================================
  ! Jacobi polynomials

  !-----------------------------------------------------------------------------
  !> Returns the Jacobi polynomial P^{a,b}_{n} at position x in [-1,1]

  pure function JacobiPolynomial(a, b, n, x) result(y)
    real(RNP), intent(in) :: a  !< first exponent, a = α
    real(RNP), intent(in) :: b  !< second exponent, b = β
    integer,   intent(in) :: n  !< order of the Jacobi polynomial
    real(RNP), intent(in) :: x  !< position in [-1,1]
    real(RNP)             :: y  !< P^{a,b}_{n}(x)
    ...
  end function JacobiPolynomial
  ...
  !=============================================================================

end module Gauss_Jacobi
```

This example demonstrates the usage of in-source documentation and illustrates further rules:

  - Apply the `only` option of the `use` statement whenever feasible.
  - Use `implicit none` to avoid errors due to implicit typing.
  - Use the `private` statement to make all module entities invisible by default
  - Use kind parameters when defining real or complex valued constants and variables, i.e.,
      * `RNP` _normal_ precision, typically 8 Byte,
      * `RHP` _high_ precision, typically 16 Byte,
      * `RDP` for interfacing libraries that require `double precision` arguments, e.g. BLAS and LAPACK.
  - Use the default kind for other intrinsic types, i.e., `integer`, `logical`, and `character`.
      

### Procedures

Functions and subroutines are preceded by a documentation header comprising a headline that is automatically adopted as the `summary` and, optionally, a detailed description. If appropriate, `date` and `author` fields  can be added as well as `todo` or other environments.

Procedure arguments are declared declared immediately after the function or subroutine statement. Further:

  - Use a separate line for each argument.
  - Specify the intent of the argument.
  - Provide a short description to each argument by adding an inline
comment.
  - If necessary, define the function result just after the arguments.
  - Insert a blank line between these declarations and the body of the procedure.
  - Put another blank line between the body and the end statement.
  - Include type and the name of the procedure in `end` statement 
  
For illustration consider the function `JacobiPolynomial` in the [`module Gauss_Jacobi`](#module_Gauss_Jacobi) and the following example of a subroutine .

<a name="subroutine_IntegerSort">_Listing:_</a> Stripped-down version of subroutine `IntegerSort` contained in `quick_sort.f`.

```Fortran  
  !-----------------------------------------------------------------------------
  !> Returns x such that x(j) <= x(k) for all j < k

  pure recursive subroutine IntegerSort(x, i1, i2)
    integer,           intent(inout) :: x(:) !< array, sorted on output
    integer, optional, intent(in)    :: i1   !< start index
    integer, optional, intent(in)    :: i2   !< terminal index

    integer :: j1, j2, j10, j20 ! internal variables
  
    ...

  end subroutine IntegerSort
```

### Types

User defined types are no genuine program units, but resemble procedures and are documented in a similar fashion. This is a simple example:

<a name="type_StandardOperators1D">_Listing:_</a> Type defining a mesh-element face.

```Fortran
  !-----------------------------------------------------------------------------
  !> Structure for keeping element face data

  type ElementFace
    integer :: id       = 0  !< mesh face ID
    integer :: neighbor = 0  !< neighbor element, including virtual one
    integer :: boundary = 0  !< adjacent boundary, if any
  end type ElementFace
```

Note that, by default, all components are `publicì. Moreover, the above type makes use of component initialization.
The next listing gives a comprehensive example.

<a name="type_StandardOperators1D">_Listing:_</a> Type providing the 1D standard-element operators, extracted from `standard_operators_1d.f`.

```Fortran
  !-----------------------------------------------------------------------------
  !> Standard operators for one polynomial order
  !>
  !> Supports the following types of nodal base functions
  !>
  !>   * Lagrange polynomials to Gauss-Legendre points:         `basis = 'GL'`
  !>   * Lagrange polynomials to Gauss-Lobatto-Legendre points: `basis = 'GLL'`
  !>   * Lagrange polynomials to Gauss-Radau-Legendre points:   `basis = 'GRL'`

  type, public :: StandardOperators1D
    private

    ! public components
    character(len=3)      , public :: basis    !< basis type
    integer               , public :: po = -1  !< polynomial order
    real(RNP), allocatable, public :: x(:)     !< collocation points
    real(RNP), allocatable, public :: w(:)     !< quadrature weights
    real(RNP), allocatable, public :: D(:,:)   !< differention matrix
    real(RNP), allocatable, public :: L(:,:)   !< stiffness (Laplace) matrix

    ! private components
    real(RNP), allocatable :: VL(:,:)     !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: VL_inv(:,:) !< inverse Legendre-Vandermonde matrix

  contains

    procedure :: Init_StandardOperators1D ! type-bound procedures
    procedure :: PolynomialOrder          ! are documented
    procedure :: InitLegendreVDM          ! separately 
    procedure :: HasLegendreVDM           ! in the same way 
    procedure :: GetLegendreVDM           ! as ordinary
    procedure :: GetInverseLegendreVDM    ! procedures

  end type StandardOperators1D
```

In this case, the default visibility is changed to `private`. Only components with the `public` attribute are accessible. These settings do not affect the type-bound procedures listed in the `contains` part. Their visibility follows identical rules, but must be declared separately.


## Programming guidelines


### Arrays and dynamic memory

Best practice rules:

  - Prefer automatic over allocatable arrays, and those over pointers.
  - In dynamic derived types prefer allocatable over pointer components.
  - Use the `associate` construct to define shorthands for deeply structured type components. 
  - Use `move_alloc` to move an allocation from one allocatable object to another.
  - Use allocatable arrays with the save attribute for reusing internal work space.

Beware that automatic arrays are often placed at the so-called stack memory, which is often limited. On Unix-like systems the stack size can be increased manually using the shell command `ulimit`, e.g.

    ulimit -s unlimited


### Types


#### Type initialization

Unlike C++ classes, Fortran types have no constructor method. Instead the type name itself gives access to an inbuilt constructor that acts like a `function` and, hence, can be used to initialize objects by assignment. The following listing provides a simple example.

<a name="type_initialization_simple_">_Listing:_</a> Initialization using the inbuilt constructor.

```Fortran
  type Point2D
    real(RNP) :: x = 0  !< x-coordinate
    real(RNP) :: y = 0  !< y-coordinate
  end type Point2D
  
  type(Point2D) :: point = Point2D(y = 1.0_RNP)
```

Here, component `x` is not specified and thus keeps its default value defined in the type declaration. Note, however, that any component with no default initialization must be defined in the constructor.

The inbuilt constructor can be supplemented or overridden by overloading user-defined functions, e.g.

<a name="constructor_overloading_">_Listing:_</a> Constructor overloading.


```Fortran
  ! in declaration part of the module, following the definition of type Point2D

  interface Point2D
    module procedure New_Point2D
  end interface 
  ...
  
contains  

  ...
  type(Point2D) function New_Point2D(input_file) result(point)
    character(len=*) :: input_file  !< input file
    
    real(RNP) :: x, y
    
    ! read x, y from input file
    ...
    
    point % x = x
    point % y = y
  
  end function New_Point2D
```

New pecularities may occur with type extension. For example consider the type `MassPoint2D` extending `Point2D`.

<a name="type_extension">_Listing:_</a> Type extension.

```Fortran
  type, extends(Point2D) :: MassPoint2D
    real(RNP) :: m  !< point mass
  end type MassPoint2D
  
  type(MassPoint2D) :: mass_point = MassPoint2D(x = 1.0_RNP, m = 10.0_RNP)
```

The child type inherits all components and, if defined, all type-bound procedures of the parent. It also possesses an inbuilt constructor similar to the parent. However, the user-defined constructor of the parent does not easily extend to the child. Instead a new, specific constructor must be provided. Nevertheless, it is possible to structure the parent constructor in such a way that it can be used to initialize the parent components within the child type. To see, how accomplish this, confer to the modules

  - `Standard_Operators_1D`, providing parent type `StandardOperators1D`, and
  - `CG_Element_Operators_1D`, defining child type `CG_ElementOperators1D`, or
  - `IP_Element_Operators_1D`, defining child type `IP_ElementOperators1D`

