\name{glinternet}
\alias{glinternet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a linear interaction model with group-lasso regularization
  that enforces 
  strong hierarchy in the estimated coefficients}
\description{The regularization path is computed along a grid of values
  for the regularization parameter lambda. Can deal with categorical
  variables with arbitrary numbers of levels, continuous variables, and
  combinations of the two. Accommodates squared error loss and logistic loss.

  The multicore option requires that the package be compiled with OpenMP support.
  Examples of compilers that qualify include gcc (>= 4.2) and icc. We also recommend
  a higher level of optimization, such as -O3 in gcc.}
\usage{
glinternet(X, Y, numLevels, lambda = NULL, nLambda = 50, lambdaMinRatio = 0.01,
 interactionCandidates=NULL, interactionPairs=NULL, screenLimit = NULL, numToFind = NULL,
family = c("gaussian","binomial"), tol = 1e-05, maxIter=5000, verbose=FALSE,
numCores = 1, weights = NULL)
}

\arguments{
  \item{X}{Matrix of features or predictors with dimension nobs x nvars;
    each row is an observation vector. Categorical variables must be coded
    as 0, 1, 2, ...}
  \item{Y}{Target variable of length nobs. Continuous for squared error
    loss, 0-1 for logistic loss.}
  \item{numLevels}{Number of levels for each variable, of length
    nvars. Set to 1 for continuous variables.}
  \item{lambda}{A user supplied \code{lambda} sequence. Typical usage
    is to have the program compute its own \code{lambda} sequence based on
    \code{nLambda} and \code{lambdaMinRatio}. Supplying a value of
    \code{lambda} overrides this. Values must be finite, positive, and strictly
    decreasing; duplicates are rejected. Every supplied value is fitted
    exactly, including when a single value is supplied.}
  \item{nLambda}{The number of \code{lambda} values. Default is 50.}
  \item{lambdaMinRatio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambdaMax}, the (data derived) entry value (i.e. the smallest
    value for which all coefficients are zero). The default is
    \code{0.01}.}
  \item{interactionCandidates}{An optional vector of variable
    indices. This will force the algorithm to only consider interactions
    between \code{interactionCandidates} and all other variables.}
  \item{interactionPairs}{An optional nx2 matrix of variable
    indices. This will force the algorithm to only consider the interaction
    pairs defined by this matrix. For example, matrix(c(1,2,1,5),
    ncol=2, byrow=TRUE) restricts the model to two interaction pairs:
    one between variables 1 and 2, and another between 1 and 5. Pair order is
    ignored and duplicate pairs are removed.}
  \item{screenLimit}{If not null (the default), limits the size of the interaction
    search space to \code{screenLimit} x nvars by only considering
    interactions with the best \code{screenLimit} candidate main effects
    at each point along the regularization path. Set this accordingly for
    large problems or if there are memory limitations. This is a heuristic
    restriction of the interaction universe; KKT checks apply to retained
    candidates, so a screened fit need not equal the unrestricted fit.}
  \item{numToFind}{Stops the program after \code{numToFind} interaction
    pairs are found. Default is null - fit all values of \code{lambda}.}
  \item{family}{A character string describing the target variable:
    "gaussian" for continuous (the default), "binomial" for logistic.}
  \item{tol}{Convergence tolerance in the adaptive FISTA algorithm.}
  \item{maxIter}{Maximum number of iterations in adaptive FISTA. Default 5000.}
  \item{verbose}{Prints progress. False by default.}
  \item{numCores}{Number of threads to run. For this to work, the
    package must be installed with OpenMP enabled. Default is 1 thread.}
  \item{weights}{Optional nonnegative case/frequency weights, one per
    observation. Weights must be finite and have positive total. They are
    normalized internally to sum to the number of observations.}
}
\details{
  The sequence of models implied by \code{lambda} is fit by FISTA (fast
  iterative soft thresholding) with adaptive step size and adaptive
  momentum restart. The continuous features are standardized to have
  weighted unit norm and weighted mean zero before computing the \code{lambda} sequence
  and fitting. The returned coefficients are unstandardized. Categorical
  variables are not standardized.

  A computed path supports \code{nLambda = 1}. If every candidate score is
  zero, the computed lambda grid is zero and the intercept-only path is
  returned. User-supplied lambda values must remain strictly positive.

  The supplied values are case/frequency weights, not analytic inverse-variance
  weights. If \eqn{w_i} is the supplied weight, the implementation uses
  \eqn{v_i = n w_i / \sum_j w_j} and minimizes
  \eqn{n^{-1}\sum_i v_i L(y_i,\eta_i) + \lambda\sum_g ||\beta_g||_2}.
  Here \eqn{L=(y-\eta)^2/2} for Gaussian data and
  \eqn{L=\log(1+\exp(\eta))-y\eta} for binomial data; the existing
  hierarchical group penalty is unchanged. Zero-weight observations are
  allowed and make no contribution to fitting, screening, KKT checks, or
  preprocessing moments. Multiplying every weight by the same positive scalar
  leaves the fitted path unchanged. Omitting weights, or supplying all-one
  weights, retains the original unweighted objective and lambda scale.
  Binomial fits require both response classes to have positive total weight.

  When \code{screenLimit} is used, it heuristically restricts the interaction
  candidate universe. Consequently, the result can differ from an unrestricted
  fit; for path entries with \code{converged = TRUE}, KKT checks and approximate
  first-order optimality apply within the retained candidates only.
}
\value{
  An object of class \code{glinternet} with the components
  \item{call}{The user function call.}
  \item{fitted}{The fitted values, with dimension nobs x
    \code{nLambda}. If \code{numToFind} is specified, the program is
    likely to stop before all \code{nLambda} models have been fit.}
  \item{lambda}{The actual \code{lambda} sequence used.}
  \item{objValue}{Objective values for each lambda.}
  \item{converged}{Logical vector indicating whether FISTA converged for each
    fitted lambda. If an iteration limit is reached, the fit is returned with
    the corresponding value set to \code{FALSE} and a warning is issued.}
  \item{iterations}{Cumulative FISTA iteration count for each fitted lambda,
    including any KKT repair solves.}
  \item{activeSet}{A list (of length \code{lambda}) of the variables found. Internally, the
    categorical and continuous variables are separated into two groups,
    and each has their own indexing system (1-based). For example, the
    categorical-continuous interaction \code{c(i, j)} refers to the
    interaction between the \code{i}-th categorical variable with the
    \code{j}-th continuous variable.}
  \item{betahat}{List (of length \code{lambda}) of coefficients for the variables in
    \code{activeSet}. The first component is the intercept. Subsequent entries correspond to the variables in activeSet. For example, if the first variable in activeSet is a 3-level categorical variable, then components 2-4 of betahat are the coefficients for this variable.}
  \item{numLevels}{The number of levels for each variable.}
  \item{family}{The target variable type.}
  \item{weights}{The internally normalized observation weights.}
}
\references{
Michael Lim and Trevor Hastie (2013)\emph{Learning interactions via hierarchical group-lasso
  regularization},
\url{https://arxiv.org/abs/1308.2719}
}
\author{Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}}

\seealso{
\code{glinternet.cv}, \code{predict.glinternet},
\code{predict.glinternet.cv}, \code{plot.glinternet.cv}, \code{coef.glinternet}
}
\examples{
# gaussian response, continuous features
Y = rnorm(100)
X = matrix(rnorm(100*10), nrow=100)
numLevels = rep(1, 10)
fit = glinternet(X, Y, numLevels)

#binary response, continuous features
Y = rbinom(100, 1, 0.5)
fit = glinternet(X, Y, numLevels, family="binomial")

#binary response, categorical variables
X = matrix(sample(0:2, 100*10, replace=TRUE), nrow=100)
numLevels = rep(3, 10)
fit = glinternet(X, Y, numLevels, family="binomial")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
