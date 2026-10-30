# User-defined INS problems

Besides the predefined problems in `src/ins_flow/problem`, the incompressible
flow solvers accept one user-defined problem per build. It is selected in the
solver input file with

```
flow_problem = 'User'
```

The user problem is compiled into the library `HiINSFlow`, hence it is
available to all programs linking this library, e.g. `ins_integrator_3d_test`
and `ml_ins_solver_3d`.

## Contents

| Path        | Purpose                            |
|-------------|------------------------------------|
| `example/`  | shipped example and template       |
| `<name>/`   | your own problems (ignored by git) |

Without user problem, the dummy `src/ins_flow/problem/ins__problem__user__3d.f`
is compiled.

## Writing a problem

A user problem is a directory with a `CMakeLists.txt` and the Fortran
sources. One source defines

- the module `INS__Problem__User__3D` with
- the type `INS_Problem_User_3D`, which extends `INS_Problem_3D`
  (`src/ins_flow/problem/ins__problem__3d.f`).

These names are fixed, since the solvers refer to them. The type must provide
the deferred procedures `SetProblem`, `GetInitialValues` and
`GetBoundaryValues`; further procedures of `INS_Problem_3D` may be overridden,
e.g. `GetExternalSources`, `GetViscosity` or the `GetExact*` family.
Problem parameters are usually read in `SetProblem` from the file given by
`problem_file` in the solver input, see the example.

Start from the example:

```
cp -r user/ins_problem/example user/ins_problem/my_flow
```

The directory may also be located outside of HiSPEET.

## CMakeLists.txt

The `CMakeLists.txt` of the problem directory calls

```cmake
hispeet_ins_user_problem(
    SOURCES   ins__problem__user__3d.f my_helpers.f
    LIBRARIES MyLib::MyLib
)
```

- `SOURCES`: all Fortran sources of the problem, relative to the directory.
  They are compiled with the same compiler options as `HiINSFlow`, including
  OpenMP and offloading, and may use all HiSPEET modules.
- `LIBRARIES` (optional): libraries used by the problem, preferably imported
  targets found with `find_package()` before the call. Plain library paths
  are accepted as well.

### Prebuilt library with module files

A Fortran library usually consists of an archive, e.g. `libmylib.a`, and the
module files (`*.mod`) of its interface. Both must be known: the archive for
linking, the module directory for compiling the sources that `use` the
library. Packages found with `find_package()` normally provide both through
an imported target. Otherwise, define the imported target yourself:

```cmake
# location of the library, e.g. -D MYLIB_DIR=/opt/mylib
set(MYLIB_DIR "" CACHE PATH "MyLib installation directory")

add_library(MyLib STATIC IMPORTED)
set_target_properties(MyLib PROPERTIES
    IMPORTED_LOCATION             "${MYLIB_DIR}/lib/libmylib.a"
    INTERFACE_INCLUDE_DIRECTORIES "${MYLIB_DIR}/include")  # *.mod files

hispeet_ins_user_problem(
    SOURCES   ins__problem__user__3d.f
    LIBRARIES MyLib
)
```

Module files are compiler specific: the library must be compiled with the
same compiler (and a compatible version) as HiSPEET. Libraries needed by
`libmylib.a` itself are added to `LIBRARIES` as well, or to the imported
target via `INTERFACE_LINK_LIBRARIES`.

## Selecting the problem

Configure HiSPEET with

```
cmake -D INS_USER_PROBLEM=<value> ...
```

| Value                | Problem directory              |
|----------------------|--------------------------------|
| empty (default)      | none, the dummy is compiled    |
| `<name>`             | `user/ins_problem/<name>`      |
| absolute path        | the given directory            |

CMake reports the choice as `INS user problem: ...`. With the dummy, selecting
`flow_problem = 'User'` stops with an error message.

To work with several problems, use one build directory per problem or
reconfigure an existing one; then only `HiINSFlow` is recompiled and the
programs are relinked.

## Example

`example/` is a copy of the predefined problem `Transition_TG` (decay of a
disturbed Taylor-Green vortex) and yields identical results. To run it, take
the input files of the predefined problem,

```
develop/ins_flow/ins_integrator_3d/transition_tg/test_mg-p.prm
develop/ins_flow/ins_integrator_3d/transition_tg/transition_tg.prm
```

and replace `flow_problem = 'Transition_TG'` by `flow_problem = 'User'`.
