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
- Fit every user-supplied lambda value, including a length-one sequence, without
  silently prepending a computed `lambdaMax`.
- Preserve matrix shapes for single-lambda fits and predictions, and return a
  well-defined intercept-only path when all candidate scores are zero.

## Quality and maintenance

- Add regression, validation, invariance, hierarchy, and cross-validation tests.
- Harden native edge cases for restricted interaction pairs, degenerate weighted
  columns, zero-weight rows, and stalled or non-finite optimization steps.
- Track convergence and iteration counts for every cross-validation fold/lambda
  fit and warn when nonconvergence can make selection unreliable.
- Preserve one-column active-set dimensions when KKT checks add multiple main
  effects, preventing violators from being recycled or omitted.
- Use scale-safe continuous standardization and reject categorical group sizes
  that exceed native integer limits before allocation.
- Compute binomial null intercepts robustly from weighted class masses, including
  extremely imbalanced positive weights.
- Replace assertion-style failures with clear validation errors for response,
  design, level, lambda, and interaction-index inputs.
- Add automated multi-platform `R CMD check` on pull requests and default-branch
  pushes.
- Document the weighted objective, input constraints, normalization, and fitted
  object components.
