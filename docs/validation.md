# Validation and simulation guide

## Automated tests

The test suite covers:

- input and weight validation;
- all-one equivalence and common-scale invariance;
- exact weighted centering/scaling geometry;
- direct reconstruction of Gaussian and binomial objectives;
- lambda-maximum behavior and fixed-lambda fitting;
- zero-weight contamination resistance;
- weighted CV aggregation and reproducible folds;
- strong hierarchy in extracted effects;
- inactive-candidate KKT checks;
- full-vector active-group KKT residuals in a tightly converged Gaussian fit;
- integer-replication equivalence with the row-geometry lambda adjustment; and
- explicit exposure of iteration-limit nonconvergence.

Run the package checks from a clean R session:

```r
devtools::test()
devtools::check()
```

## Target-risk simulation

`inst/validation/weighted-target-simulation.R` is a deterministic simulation
showing why case weights matter. A target population contains two latent strata
with different treatment-modifier slopes. Sampling deliberately overrepresents
the high-modification stratum. Because stratum is omitted from the fitted model,
the interaction coefficient is a target-distribution average:

\[
  \tau_{target}=E_{target}(1+2S)=1.5.
\]

The unweighted sample targets the distorted sample mixture, whereas inverse-
sampling weights recover the target mixture. Each replicate fits the same
hierarchical treatment-by-continuous-modifier model with and without weights and
extracts its interaction slope using prediction contrasts. The script reports
bias, root mean squared error, and the proportion of replicates closer to the
target under weighting.

With the fixed seed and 30 replicates in the script, the mean interaction was
2.231 unweighted and 1.492 weighted, versus the target 1.500. Corresponding RMSE
values were 0.737 and 0.097, and the weighted estimate was closer in every
replicate.

Run it after installing the package:

```r
source(system.file(
  "validation",
  "weighted-target-simulation.R",
  package = "glinternet"
))
```

This simulation demonstrates correct influence on the fitted target risk. It is
not a universal performance guarantee: highly variable or misspecified weights
can increase variance, and causal interpretation still requires exchangeability,
positivity, consistency, and correct definition of the target population.

## Numerical certification

The automated full-vector KKT test uses a tight optimizer tolerance and a typical
well-conditioned design. It guards against broad regressions but does not erase
the source-level convergence limitation described in the audit. For critical
fits, reconstruct the native design, calculate the maximum active-group vector
residual and inactive-group norm violation, and retain the diagnostic with the
analysis artifacts.
