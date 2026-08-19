\name{glinternet.cv}
\alias{glinternet.cv}
\title{Cross-validation for glinternet}
\description{Does k-fold cross validation for glinternet and returns a
  value of \code{lambda}.}
\usage{
glinternet.cv(X, Y, numLevels, nFolds = 10, lambda=NULL, nLambda=50,
lambdaMinRatio=0.01, interactionCandidates=NULL, interactionPairs=NULL,
screenLimit=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000,
verbose=FALSE, numCores=1, weights=NULL, foldid=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\code{X} matrix as in \code{glinternet}.}
  \item{Y}{Target \code{Y} as in \code{glinternet}.}
  \item{numLevels}{Number of levels \code{numLevels} as in \code{glinternet}.}
  \item{nFolds}{Number of folds - default is 10.}
  \item{lambda}{\code{lambda} as in \code{glinternet}.}
  \item{nLambda}{\code{nLambda} as in \code{glinternet}.}
  \item{lambdaMinRatio}{\code{lambdaMinRatio} as in \code{glinternet}.}
  \item{interactionCandidates}{\code{interactionCandidates} as in
    \code{glinternet}.}
  \item{interactionPairs}{\code{interactionPairs} as in \code{glinternet}.}
  \item{screenLimit}{\code{screenLimit} as in \code{glinternet}.}
  \item{family}{\code{family} as in \code{glinternet}.}
  \item{tol}{\code{tol} as in \code{glinternet}.}
  \item{maxIter}{\code{maxIter} as in \code{glinternet}.}
  \item{verbose}{\code{verbose} as in \code{glinternet}.}
  \item{numCores}{\code{numCores} as in \code{glinternet}.}
  \item{weights}{Optional nonnegative case/frequency weights as in
    \code{glinternet}.}
  \item{foldid}{Optional vector of integer fold labels, one per observation.
    At least two nonempty folds are required and each fold must have positive
    total validation weight.}
}
\details{The \code{lambda} sequence is computed using all the
  data. \code{nFolds} models are fit, each time with one of the folds
  omitted. Each training fit renormalizes its own training weights. Validation
  losses are weighted by observation weight and combined by fold weight mass,
  so \code{cvErr} is the observation-weighted loss across all out-of-fold
  predictions. A weighted standard error across fold means is also computed;
  with equal fold masses this reduces to the fold-loss standard deviation
  divided by the square root of the number of folds. The one-standard-error
  choice is the largest lambda with error no greater than the minimum error plus
  the standard error at the minimum. The full lambda path is fit to the full
  data set, and the selected components are extracted at \code{lambdaHat}.
  If \code{screenLimit} is set, its heuristic retained-candidate
  semantics are the same as in \code{glinternet}. For binomial models, use
  suitably stratified folds so both response classes have positive total
  training weight in every fold.}
\value{
  An object of class \code{glinternet.cv} with the components
  \item{call}{The user function call.}
  \item{glinternetFit}{Glinternet object containing the full lambda path fitted
    to all observations.}
  \item{fitted}{Vector for fitted values (same length as \code{Y}). This
  is from the model fitted at \code{lambdaHat}.}
  \item{activeSet}{\code{activeSet} is a list variables found for the
    model fitted with \code{lambdaHat}.}
  \item{betahat}{Unstandardized coefficients for the variables in
    \code{activeSet}.}
  \item{lambda}{The actual sequence of \code{lambda} values used for the
    cross validation.}
  \item{lambdaHat}{The value of \code{lambda} that minimizes the cv
    error curve.}
  \item{lambdaHat1Std}{The largest value of \code{lambda} that produces
    a cv error that is within 1 standard error of the minimum cv
    error. This will always be at least as large as \code{lambdaHat}.}
  \item{cvErr}{The vector of cross validation errors. Same length as
    \code{lambda}.}
  \item{cvErrStd}{Weighted standard error for cross-validation errors across
    the \code{nFolds} fold means.}
  \item{family}{The response type.}
  \item{numLevels}{Input number of levels for each variable.}
  \item{nFolds}{The number of folds used.}
  \item{foldid}{The fold assignment used.}
  \item{foldLoss}{Matrix of per-fold mean validation losses.}
  \item{foldWeight}{Total raw validation weight in each fold.}
  \item{weights}{The raw observation weights supplied to cross-validation.}
}

\author{
Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}
}

\seealso{
\code{glinternet}, \code{predict.glinternet},
\code{predict.glinternet.cv}, \code{plot.glinternet.cv}
}
\examples{
Y = rnorm(100)
numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1)
rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
fit = glinternet.cv(X, Y, numLevels, nFolds=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
