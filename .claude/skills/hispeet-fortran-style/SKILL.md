---
name: hispeet-fortran-style
description: Fortran formatting and naming rules of HiSPEET (Jörg Stiller). Apply to every new or modified Fortran source, include file or test program (.f, .F, .f90, .F90, .inc, .F.in) before showing it for review or committing.
---

# HiSPEET Fortran style

Apply these rules to all code you write or change, not only to the lines of
the current task. Check the touched blocks once more before handing files to
Jörg for review.

## Continuation lines

- Vertically align the continuation signs `&`.
- The closest `&` is placed two spaces away from the text (the longest
  statement line of the block).
- Related blocks (e.g. consecutive calls with the same argument structure,
  OpenMP directive continuations of one construct) share one `&` column.
- Argument lists:
  - first line: one space between the opening `(` and the first argument
  - following lines: the leading `,` is aligned with the `(`
  - last line: the closing `)` is aligned with the `&` above

  ```fortran
  call VolumeKernel_R( eop % po, mesh % n_elem_active, size(fw,4)  &
                     , nm, .true., gf, gd, eop % w, eop % D        &
                     , metrics % a, metrics % Ji_n, fw, r, f       )
  ```

## OpenMP directive continuations

- The same rules apply to continued directives (`!$omp ... &`): all `&` of
  one directive in one column, two spaces after the longest continued line
  of that directive; clauses on the following lines indented by two extra
  spaces after `!$omp`.

  ```fortran
  !$omp target teams distribute parallel do collapse(4)  &
  !$omp   private(s, mf, a_0, a_P, c_0, c_P)             &
  !$omp   map(to: gf, gd, a, Ji_n)                       &
  !$omp   map(tofrom: r(:,:,:,1:ne)) map(to: f(:,:,:,1:ne))
  ```

## Columns in lists of similar statements

- In runs of similar statements (declarations, assignments, calls such as
  `call XMPI_Bcast(x, 0, comm)`), align `=`, `::`, trailing `!` comments and
  the commas after the varying first argument in one column, padded to the
  longest entry.
- Array dimensions are normally given after the variable name. Only where
  this clearly hurts readability, e.g. a long array spec in an aligned
  declaration run, a `dimension(...)` attribute may be used instead, so the
  names and `!<` comments stay in their columns:

  ```fortran
  integer,                         intent(in) :: method  !< solver {1,2,3}
  real(RNP), dimension(:,:,:,:,:), intent(in) :: x       !< mesh points
  ```

## Module layout

- Procedures after `contains` of a module are indented by two columns,
  including their header comments and the separator lines.
- Separator lines (`!---...`, `!===...`) end at column 80, also when
  indented (`  !----...` is two characters shorter than at column 1).

## Names

- Variable names separate a trailing index or qualifier by an underscore:
  `time_0`, `a_0`, `r_l2_0` (not `time0`).
- Procedures that initiate nonblocking communication are named `Start...`
  (`StartReceives`, `StartSends`), not `Post...`.

## General

- Lines at most ~80 characters.
- Fixed parameters (named constants) in UPPERCASE with descriptive names.
- Avoid code duplication: element-bound routines build on flat counterparts
  or share code via include files.
- Code, comments and commit messages in English.
