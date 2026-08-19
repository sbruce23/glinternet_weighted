# glinternet 1.0.14

## Observation weights

- Add nonnegative observation weights to Gaussian and binomial hierarchical
  group-lasso fits.
- Apply weights consistently to preprocessing, lambda-path construction,
  screening, optimization, KKT checks, coefficient rescaling, and prediction.
- Add weighted cross-validation and reproducible user-supplied fold assignments.
- Preserve unweighted behavior when weights are omitted or all equal.
- Report per-lambda convergence and iteration metadata, and warn when FISTA
  reaches `maxIter` before convergence.

## Quality and maintenance

- Add regression, validation, invariance, hierarchy, and cross-validation tests.
- Harden native edge cases for restricted interaction pairs, degenerate weighted
  columns, zero-weight rows, and stalled or non-finite optimization steps.
- Add automated multi-platform `R CMD check` on pull requests and default-branch
  pushes.
- Document the weighted objective, input constraints, normalization, and fitted
  object components.
