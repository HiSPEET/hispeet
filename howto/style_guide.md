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

<a name="module_Gauss_Jacobi">_Code example:_</a> Stripped-down version of the `Gauss_Jacobi` module contained in `gauss_jacobi.f`.

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

<a name="subroutine_IntegerSort">_Code example:_</a> Stripped-down version of subroutine `IntegerSort` contained in `quick_sort.f`.

```Fortran  

!-------------------------------------------------------------------------------
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

User defined types are no genuine program units, but in many ways similar and sometimes even more complex. In the simplest case, the type just bundles a number of components, e.g.

```Fortran
type Foo
  real :: a
end type Foo
```


