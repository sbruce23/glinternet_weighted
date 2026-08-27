# Mathematical and statistical audit

## Executive assessment

The weighted extension has a coherent convex objective and applies normalized
observation weights consistently to the central fitting operations. The existing
foundation is suitable for weighted exploratory hierarchical interaction
selection, subject to the qualifications below. The audit found one numerical
certification defect that should be fixed in a later source-code change: the
native active-group convergence check verifies gradient *norm* but not gradient
*direction*. This repository update deliberately documents and tests around that
finding without changing `R/` or `src/`.

## Verified properties

The following properties were checked by direct source inspection and independent
tests.

- Weights are finite, nonnegative, have positive total, and are normalized to
  sum to the retained row count.
- Gaussian and binomial losses, intercept updates, weighted moments,
  standardization, lambda scores, screening gradients, FISTA gradients,
  objective values, inactive-group KKT checks, and validation loss all use the
  normalized weights.
- All-one weights reproduce the legacy unweighted path to numerical precision;
  multiplying all weights by a positive constant leaves the path unchanged.
- Zero-weight values are excluded from weighted arithmetic, including safeguards
  for otherwise contaminating extreme values.
- Fixed-lambda objectives can be independently reconstructed from returned
  residuals and native group coefficients.
- The Lim–Hastie overlapping groups retain strong hierarchy in the extracted
  effects.
- Weighted cross-validation aggregates fold loss by validation weight mass and
  exposes per-fold nonconvergence.
- Integer replication agrees after the row-geometry lambda adjustment
  \(\lambda\sqrt{n/N}\).

These checks support the mathematical soundness of the weighted objective and
the claim that observation weights genuinely influence estimation.

## Finding 1: active-group convergence is not a full KKT certificate

For a nonzero group \(g\), first-order optimality requires the vector equation

\[
  \nabla_g\ell(\widehat\beta)
  +\lambda\frac{\widehat\beta_g}
                  {\lVert\widehat\beta_g\rVert_2}=0.
\]

For a zero group it requires
\(\lVert\nabla_g\ell(\widehat\beta)\rVert_2\leq\lambda\).
The native function `check_convergence()` currently checks only

\[
  \left|\lVert\nabla_g\ell\rVert_2-\lambda\right|/\lambda
  \leq \texttt{tol}
\]

for an active group. Equal norms are necessary but do not establish that the two
vectors point in opposite directions. Stress tests can therefore produce
`converged = TRUE` while the full active-group vector residual is materially
larger than the requested tolerance.

### Consequence

This is a flaw in the meaning of the convergence flag, not evidence that every
fit is wrong. Typical fits with a tight tolerance passed independent full-vector
checks. Nevertheless, a publication-quality numerical claim should not treat
`converged = TRUE` as sufficient KKT certification until the source criterion is
replaced by a proximal-gradient or full-vector KKT residual.

### Recommended later source fix

For each active group, evaluate the vector residual above, scaled by a stable
denominator such as \(\max(1,\lambda,\lVert\nabla_g\ell\rVert_2)\). For zero
groups, retain the norm inequality. Include intercept stationarity and report the
maximum residual. That source change is outside the scope of this documentation-
only audit branch.

## Finding 2: “frequency weights” needs a lambda qualification

The weighted loss is mathematically well defined. However, literal replication
changes native group scaling through row count. Same-lambda replication is not
an invariance of this package. The exact comparison uses
\(\lambda_{rep}=\lambda\sqrt{n/N}\), as documented in
[weighted-methodology.md](weighted-methodology.md#case-importance-weights-and-frequency-replication).
The accurate general description is *case/importance weights*.

## Finding 3: cross-validation uncertainty is heuristic

The CV point estimate is the correct weighted mean of out-of-fold loss. The
reported `cvErrStd` summarizes between-fold dispersion using effective fold
count. Fold losses overlap in their training data and are not independent sample
means, so this quantity should be called a fold-dispersion standard-error
heuristic rather than a formal standard error for generalization risk. This is
common and useful for a one-standard-error selection rule, but its inferential
meaning should not be overstated.

## Finding 4: screening changes the interaction universe

`screenLimit` is a computational heuristic. The inactive KKT calculation applies
to retained candidates, not necessarily to every omitted interaction. A screened
fit and an unrestricted fit need not coincide. Analyses should report the
candidate universe, screening settings, and whether conclusions persist without
screening when that comparison is computationally feasible.

## Finding 5: estimated-weight uncertainty is external

The package correctly treats supplied weights as fixed inputs. Inverse-
probability weights estimated from the same data add uncertainty and depend on
positivity, correct propensity specification, and the chosen target estimand.
Those design assumptions cannot be verified by an optimizer. Applied work should
report weight construction, balance, overlap, range, truncation, and effective
sample size; inference should use a method that accounts for estimated weights
when warranted.

## Defensible-use checklist

Before interpreting a weighted fit:

1. State the target population and why the supplied weights identify it.
2. Report raw and normalized weight distributions, arm-specific effective sample
   sizes, overlap, and post-weighting covariate balance.
3. Freeze the eligible interaction universe and disclose any screening.
4. Supply treatment-stratified, weight-balanced fold IDs for observational
   treatment-modifier work.
5. Compare minimum-CV and one-standard-error selections and assess fold/seed
   stability.
6. Require all convergence flags to be true and independently check full-vector
   KKT residuals for reported path points.
7. Treat selected interactions as exploratory unless evaluated honestly in new
   data or by a prespecified sample-splitting/resampling design.
8. Avoid interpreting penalized coefficient magnitude as an unbiased effect
   estimate; estimate clinically meaningful contrasts separately.

## Overall conclusion

The weighted formulation is built on a solid statistical foundation and is not
invalidated by the audit. The major qualification is narrower but important:
the current native convergence flag is not a mathematically complete certificate.
Documentation, independent KKT validation, tight tolerances, and sensitivity
checks make current exploratory use defensible; a later source patch is needed
before claiming that every `converged = TRUE` solution has been fully KKT-
certified.
