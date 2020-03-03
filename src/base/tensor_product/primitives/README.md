
## Notation

### Operators and operands

   Symbol  |  Usage
  -------- | ---------------------------
   `I`     |  1D unit matrix
   `A,B,C` |  1D rectangular matrices
   `D`     |  1D diagonal matrix
   `Q`     |  1D square matrix
   `S`     |  1D symmetric matrix
   `W`     |  3D diagonal matrix
   `t`     |  transposed, i.e. `Qt = Qᵀ`

### Loop transformations

Example `IxIxQt_01_01_41_44`

 * routine for evaluating `v = α I⊗I⊗Qᵀ u + β v`
 * number blocks indicate blocking and unrolling of `k`, `j`, `i` and `p` (reduction) loops
 * here
   - `k`: blocking `0`, unroll width `1` (i.e. no blocking and no unrolling)
   - `j`: blocking `0`, unroll width `1`
   - `i`: blocking `4`, unroll width `1`
   - `p`: blocking `4`, unroll width `4`

### Classification of complex operators

Names `TPO_<Name>_[R|D][L][C|V][I|A][_<dim_u><dim_v>]`

   Specifier |  Meaning
  ---------- | -------------------------------------
   `R`/`D`   |  regular/deformed elements
   `L`       |  diagonal (lumped) mass matrix
   `C`/`V`   |  constant/variable coefficients
   `I`/`A`   |  isotropic/anisotropic coefficients

Optionally, the dimensions of the operand, `dim_u`, and of the result, `dim_v`, can be specified as follows

   Specifier  |  Meaning
  ----------- | -----------------
   `S`        |  scalar
   `V`        |  vector

Typically, the specifiers will be hidden in favor of generic interfaces.