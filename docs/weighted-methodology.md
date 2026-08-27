# Weighted hierarchical group LASSO methodology

## Scope and estimand

For observations \((y_i,x_i,w_i)\), `glinternet` estimates the minimizer of a
weighted empirical risk over a strongly hierarchical pairwise-interaction model.
The nonnegative weights describe the target empirical distribution: a row with
twice the weight contributes twice as much loss. They may be survey, sampling,
inverse-probability, transport, or other importance weights when the scientific
design justifies that interpretation. The optimizer itself does not estimate a
propensity model or account for uncertainty from estimated weights.

Define normalized weights

\[
  v_i = \frac{n w_i}{\sum_{j=1}^n w_j},
  \qquad \sum_{i=1}^n v_i=n.
\]

This normalization makes a common positive rescaling of all supplied weights
irrelevant and retains the historical `glinternet` lambda scale for all-one
weights.

## Model and objective

Let \(\eta_i=\beta_0+f_\beta(x_i)\), where \(f_\beta\) contains main effects and
eligible pairwise interactions in the overlapping group parameterization of Lim
and Hastie. For Gaussian outcomes, the fitted path minimizes

\[
  Q_\lambda(\beta_0,\beta)
  = \frac{1}{2n}\sum_{i=1}^n v_i(y_i-\eta_i)^2
    +\lambda\sum_{g\in\mathcal G}\lVert\beta_g\rVert_2.
\]

For binomial outcomes, the loss term is

\[
  \frac{1}{n}\sum_{i=1}^n v_i
  \left\{\log(1+\exp(\eta_i))-y_i\eta_i\right\},
\]

implemented with a numerically stable softplus calculation. The intercept is
unpenalized. The group penalty is convex in the expanded parameterization.

## Strong hierarchy

The overlap construction duplicates lower-order coefficients inside interaction
groups. After coefficients are mapped back to the original scale, a nonzero
interaction is accompanied by both associated main effects. This is *strong
hierarchy*. It is a structural property of the parameterization, not a claim
that each selected coefficient is statistically significant.

Categorical main effects use one coefficient per level. A continuous-continuous
group contains the two parent main-effect copies and their product; mixed and
categorical-categorical groups analogously contain parent and interaction
components. Group-specific scaling follows the original package so that groups
with different dimensions are put on a comparable penalty scale.

## Weighted preprocessing

Continuous columns are centered and scaled using \(v_i\). A standardized column
has weighted mean zero and weighted squared norm one under the package's native
geometry. Continuous products are likewise weighted-centered and scaled.
Categorical encodings retain the original row-count scaling. Zero-weight rows do
not contribute through their values, but the retained row count still affects
group geometry.

The maximum lambda, sequential strong-rule scores, inactive-group KKT screen,
smooth-loss gradient, intercept update, and objective all use the normalized
weights. Consequently, weighting is integrated throughout the fit rather than
applied only after coefficient estimation.

## Case-importance weights and frequency replication

The implementation is naturally described in terms of case or importance
weights. Suppose integer weights are expanded into \(N=\sum_i w_i\) literal
rows. Replication changes the row-count factors in the categorical and
interaction design blocks. With the package's group geometry, the corresponding
replicated-data penalty is

\[
  \lambda_{\mathrm{rep}}
  = \lambda_{\mathrm{weighted}}\sqrt{n/N}.
\]

Thus integer weights reproduce literal replication after this lambda adjustment,
not at the same numeric lambda. This distinction does not undermine the weighted
empirical-risk formulation, but it means documentation should not promise
same-lambda frequency-weight equivalence.

## Weighted cross-validation

Within training fold \(-k\), weights are renormalized to the training row count
and the complete lambda sequence is fitted. Let \(W_k=\sum_{i\in k}w_i\) and
let \(L_k(\lambda)\) be validation-fold mean loss normalized by \(W_k\). The
reported cross-validation curve is

\[
  \bar L(\lambda)=\sum_{k=1}^K a_k L_k(\lambda),
  \qquad a_k=\frac{W_k}{\sum_j W_j}.
\]

This equals the weighted mean loss across all out-of-fold predictions. The
minimum-CV rule selects the lambda with the smallest \(\bar L\). The
one-standard-error rule uses the largest eligible lambda.

The reported uncertainty curve is a weighted dispersion of the fold means:

\[
  K_{\mathrm{eff}}=\frac{1}{\sum_k a_k^2},\qquad
  s^2(\lambda)=
  \frac{\sum_k a_k\{L_k(\lambda)-\bar L(\lambda)\}^2}
       {1-\sum_k a_k^2},\qquad
  \mathrm{SE}_{\mathrm{fold}}=\sqrt{s^2/K_{\mathrm{eff}}}.
\]

It is a useful fold-dispersion heuristic for model selection, not a formal
sampling standard error for prediction risk. When weights are estimated, neither
this quantity nor the fitted object includes propensity-model uncertainty.

## Statistical interpretation

The result is a regularized prediction and interaction-selection procedure.
Ordinary coefficient p-values are not produced, and fitting an unpenalized model
to selected variables in the same sample does not create valid post-selection
inference. Scientific applications should distinguish:

1. the target population encoded by the weights;
2. exploratory interaction discovery by penalized weighted risk;
3. validation or inference using an independent sample, honest sample split, or
   a resampling procedure that repeats every data-adaptive step; and
4. uncertainty from estimating the weights, preferably handled by stacked
   estimating equations or a bootstrap that refits the weight model.

## Reference

Michael Lim and Trevor Hastie (2015), “Learning interactions via hierarchical
group-lasso regularization,” *Journal of Computational and Graphical Statistics*,
24(3), 627–654. <https://doi.org/10.1080/10618600.2014.938812>
