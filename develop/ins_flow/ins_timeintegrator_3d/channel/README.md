# Turbulent channel flow

## TBD

- add description

## Tips for increasing robustness

- Activate filtering of intermediate velocity after extrapolation step (accuracy?)
- Use coupled projection-diffusion solver (cost?)
- Increase  `mu_0` and/or  `i_max_v`, `i_pre_v`, `i_krylov` (cost?)
- Use equal order for velocity and pressure (inf-sup stability?)

Note that every option may have adverse side effects, which are indicated in brackets. 