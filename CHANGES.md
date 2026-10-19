# Changes in branch `hpc-course`

Base: `master` at `82d3550d` (2026-10-15, "LIBXSMM: update to c14cbc6; link libxsmmext only if built").

This branch adds the exercise templates of the course *Simulationstechnik in der Strömungsmechanik* (SEM and HPC parts) to the current HiSPEET code. The templates come from branch `hpc-course` (`94161b1a`, 2022-12-06) of
`gitlab.hrz.tu-chemnitz.de/sim-tech-stroem/hispeet_exercise_solutions`.

Only `app/exercises/` differs from `master`. The libraries in `src/`, `develop/`, `external/` (including the LIBXSMM integration), the CMake setup, the README and the wiki are identical to `master`.

## Exercises

All exercise programs are declared `EXCLUDE_FROM_ALL`. A plain `make` builds the HiSPEET libraries but no exercises; build them by name, e.g. `make sem_basics`. They link against `HiBase` only.

| Folder | Exercise sheet | Programs | Other files |
|---|---|---|---|
| `ex_1` | Intro | `hello`, `calculate_pi`, `jacobi_method` (new) | `jacobi_method.py` |
| `ex_2` | SEM 1 | `sem_basics` (was `exercise_1`) | `README.md`, `sem_basics.prm`, `sem_basics.py` |
| `ex_3` | SEM 2 | `helmholtz_1d` | `README.md`, `helmholtz_1d.prm`, `helmholtz_1d.py`, `helmholtz_1d_study.py` |
| `ex_4` | SEM 3 | `exact_integration`, `tensor_products` | `README.md`, `tensor_products.prm`, `exact_integration.py`, `tensor_products.py` |
| `ex_5` | HPC 1 | `stencil_test`, `matrix_matrix_test`, `potpourri_1`, `potpourri_2`, `potpourri_3` | `README.md`, `potpourri_print_mlups.f` |
| `ex_6` | HPC 2 | `test_interpolation`, `potential_flow_solver` | `spectral_element_interpolation.f`, `helmholtz_2d.f`, input files |
| `ex_7` | HPC 3 | `sem_1d`, `ping_pong` | `input_sem_1d` |

HPC 4 is a pen-and-paper exercise and has no code.

## File-level changes compared to `master`

- `app/exercises/CMakeLists.txt`: also adds `ex_3` to `ex_7` (`master` has only `ex_1` and `ex_2`).
- `ex_1/CMakeLists.txt`: adds the `jacobi_method` target.
- `ex_2/CMakeLists.txt`: builds `sem_basics` instead of `exercise_1`.
- `ex_2/exercise_1.f` renamed to `ex_2/sem_basics.f` and reworked.
- New: all other files listed in the table above (41 files changed in total, +3550/−81 lines).

## Adaptations to the current HiSPEET

- `ex_1/CMakeLists.txt`: `hello` is linked to `HiBase` (`target_link_libraries(hello HiBase)`), so the intro task (`use Constants`, `PI`, `real(RNP)`) compiles. Before, `make hello` failed with "Cannot open module file 'kind_parameters.mod'".
- `ex_3/README.md`: the study plot is opened as `helmholtz_1D_study.pdf` (capital D), the name `helmholtz_1d_study.py` actually writes. The exercise sheet `ex_sem2.tex` still says `helmholtz_1d_study.pdf`.
- `ex_6/helmholtz_2d.f`, `ex_7/sem_1d.f`: `GLL_Points`, `GLL_Weights`, `GLL_DiffMatrix` replaced by `LobattoPoints`, `LobattoWeights`, `LobattoDiffMatrix` (same arguments). These two files did not compile even on the original 2022 branch.
- `ex_4/tensor_products.py`: `fig.savefig(outfile, layout='tight')` replaced by `fig.tight_layout()` + `fig.savefig(outfile)` (the `layout` argument raises a `TypeError` from Matplotlib 3.6 on, i.e. with the packages of Ubuntu 24.04 to 26.04; older versions silently ignored it); `usetex=True` removed, so no LaTeX installation is needed (labels use Matplotlib's built-in math rendering); non-interactive `Agg` backend set, so the script also works without a display. Tested with Matplotlib 3.1.3, 3.5.1, 3.6.3 and 3.10.8.

## Merge notes

The merge of the 2022 branch into the current code had conflicts in `CMakeLists.txt`, `src/base/CMakeLists.txt`, `external/cmake/libxsmm.cmake`, `README.md`, `src/base/time_integration/imex_runge_kutta_method.f` and `howto/wiki/ubuntu.md`. All were resolved in favour of `master`. The 2022 branch also commented out three post-processing programs in `develop/cart/isp_flow/` (`channel_post`, `change_polynomial_order`, `transition_tg_post`); they stay enabled.

## Verification

All 16 programs build with GCC 13 / Open MPI and produce the same results as on the original 2022 branch. Output files are identical, except for timings and round-off at the 1e-16 level in `tensor_products.dat`.

## Known issues (not fixed yet)

- `ex_5`: the results of the timed loops are never used, so at the default `-O3` the compiler removes them and the reported times/MLUPS are meaningless.
- `-D OPT=...` has no effect with GNU compilers (`CMakeLists.txt`, GNU section), so the optimisation-level comparison in HPC 1 cannot be done as described.
- OpenMP is off by default; HPC 2 and the OpenMP part of HPC 3 need `cmake -D OpenMP=1 ..`. Without it the `!$omp` directives are silently ignored and there is no speed-up. **To do:** add this configure step to the exercise sheets `ex_hpc2` (and its slides) and `ex_hpc3`.
