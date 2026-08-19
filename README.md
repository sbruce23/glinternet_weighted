# glinternet

`glinternet` fits pairwise interaction models with hierarchical group-lasso
regularization. It supports Gaussian and binomial responses, continuous and
categorical predictors, and strong hierarchy: an estimated interaction is
accompanied by both associated main effects.

This version also supports an observation weight for every row. Omit `weights`
for an ordinary unweighted fit.

## Installation

Install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("sbruce23/glinternet_weighted")
```

The package contains compiled C code. A suitable compiler toolchain is required
when installing from source. OpenMP is optional and enables `numCores > 1` on
platforms that support it.

## Weighted fit

```r
library(glinternet)

set.seed(1)
n <- 200
X <- matrix(rnorm(n * 6), nrow = n)
Y <- 1 + X[, 1] - X[, 2] + 2 * X[, 1] * X[, 2] + rnorm(n)
w <- runif(n, 0.5, 2)

fit <- glinternet(
  X,
  Y,
  numLevels = rep(1, ncol(X)),
  weights = w,
  family = "gaussian"
)
```

Weights are case/frequency weights: they control each row's contribution to the
likelihood. They are not inverse-variance analytic weights. Weights must be
finite and nonnegative and must have a positive total. Zero weights are allowed
and those rows make no contribution to fitting, screening, preprocessing
moments, KKT checks, or validation loss. Internally, weights are normalized to
sum to the number of observations, so multiplying all supplied weights by the
same positive constant does not change the fitted path.

For normalized weights \(\widetilde w_i = n w_i / \sum_i w_i\), the fitted
Gaussian path minimizes

\[
  \frac{1}{2n}\sum_{i=1}^n \widetilde w_i
  (y_i - \eta_i)^2 + \lambda \sum_g \lVert \beta_g \rVert_2,
\]

and the binomial path replaces squared error with binomial negative
log-likelihood,
\(L(y_i,\eta_i)=\log(1+\exp(\eta_i))-y_i\eta_i\), evaluated numerically
stably. The groups and overlapping parameterization are those of Lim and Hastie
and enforce strong hierarchy. Continuous columns and continuous-continuous
products are centered and scaled using these normalized weights.

## Weighted cross-validation

```r
foldid <- sample(rep(seq_len(5), length.out = n))
cvfit <- glinternet.cv(
  X,
  Y,
  numLevels = rep(1, ncol(X)),
  weights = w,
  nFolds = 5,
  foldid = foldid,
  family = "gaussian"
)

cvfit$lambdaHat
cvfit$lambdaHat1Std
```

Each training fold normalizes its own training weights. Within validation fold
\(k\), the mean loss is weighted by the raw validation weights. Fold losses are
then combined in proportion to their validation weight masses, which is the
weighted loss across all out-of-fold predictions. Every validation fold must
have positive weight mass.

`cvErrStd` is a standard error across fold means. If
\(a_k=W_k/\sum_j W_j\), the implementation uses effective fold count
\(K_{eff}=1/\sum_k a_k^2\), weighted variance
\(s^2=\sum_k a_k(L_k-\bar L)^2/(1-\sum_k a_k^2)\), and
\(\operatorname{SE}=\sqrt{s^2/K_{eff}}\). With equal fold masses this is
`sd(fold_loss) / sqrt(nFolds)`. The one-standard-error choice is the largest
lambda whose loss is no more than the minimum loss plus the standard error at
the minimum. Supplying `foldid` makes the split reproducible and enables direct
comparisons across fits.

For binomial cross-validation, choose stratified folds so that both response
classes have positive total training weight in every fold.

`weights = NULL` and all-one weights retain the original unweighted objective
and lambda scale. Rescaling all weights by the same positive constant also
leaves the result unchanged.

`screenLimit` retains the package's original heuristic behavior: it restricts
the interaction universe for speed and memory use. KKT checks cover retained
candidates, so screened and unrestricted fits are not guaranteed to match.

See `?glinternet` and `?glinternet.cv` for the full API and returned objects.

## Reference

Michael Lim and Trevor Hastie (2015), “Learning interactions via hierarchical
group-lasso regularization,” *Journal of Computational and Graphical
Statistics*, 24(3), 627–654. <https://doi.org/10.1080/10618600.2014.938812>
