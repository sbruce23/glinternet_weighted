# glinternet 1.0.14

## Observation weights

- Add nonnegative observation weights to Gaussian and binomial hierarchical
  group-lasso fits.
- Apply weights consistently to preprocessing, lambda-path construction,
  screening, optimization, KKT checks, coefficient rescaling, and prediction.
- Add weighted cross-validation and reproducible user-supplied fold assignments.
- Preserve unweighted behavior when weights are omitted or all equal.

## Quality and maintenance

- Add regression, validation, invariance, hierarchy, and cross-validation tests.
- Add automated multi-platform `R CMD check` on pull requests and default-branch
  pushes.
- Document the weighted objective, input constraints, normalization, and fitted
  object components.
