test_that("all five weighted candidate scores match direct design gradients", {
  set.seed(921)
  n <- 30
  Xcat <- cbind(sample(0:1,n,TRUE), sample(0:2,n,TRUE))
  Xraw <- cbind(rnorm(n), rnorm(n))
  v <- glinternet:::validate_weights(runif(n,.2,2), n)
  Z <- apply(Xraw, 2, glinternet:::standardize, weights=v)
  r <- rnorm(n)
  cand <- glinternet:::get_candidates(Xcat, Z, r, v, n, 2, 2, c(2,3),
                                      NULL, NULL, NULL, NULL)
  wr <- v*r
  dummy <- function(x, L) outer(x, 0:(L-1), `==`) + 0
  norm_score <- function(D) sqrt(sum(drop(crossprod(D, wr)/n)^2))

  expect_equal(cand$norms$cat[1], norm_score(dummy(Xcat[,1],2)/sqrt(n)), tolerance=1e-12)
  expect_equal(cand$norms$cont[1], norm_score(Z[,1,drop=FALSE]), tolerance=1e-12)
  joint <- dummy(Xcat[,1]+2*Xcat[,2], 6)/sqrt(n)
  expect_equal(cand$norms$catcat[1], norm_score(joint), tolerance=1e-12)

  u <- Z[,1]*Z[,2]
  u <- u-sum(v*u)/n
  unorm <- sqrt(sum(v*u^2))
  us <- if (unorm > 1e-15) u/unorm else rep(0,n)
  Dcc <- cbind(Z[,1],Z[,2],us)/sqrt(3)
  expect_equal(cand$norms$contcont[1], norm_score(Dcc), tolerance=1e-11)

  Dmix <- cbind(dummy(Xcat[,1],2)/sqrt(2*n),
                dummy(Xcat[,1],2)*Z[,1]/sqrt(2))
  expect_equal(cand$norms$catcont[1], norm_score(Dmix), tolerance=1e-11)
})

test_that("lambda maximum is intercept-only and solver objective is independent", {
  set.seed(193)
  n <- 45
  X <- cbind(sample(0:1,n,TRUE), rnorm(n), rnorm(n))
  y <- .4+X[,2]+rnorm(n)
  v <- glinternet:::validate_weights(runif(n,.2,2), n)
  fit <- glinternet(X,y,c(2,1,1),nLambda=4,weights=v,maxIter=4000)
  expect_true(all(vapply(fit$activeSet[[1]], is.null, logical(1))))

  Xcat <- matrix(X[,1], ncol=1)
  Z <- apply(X[,2:3,drop=FALSE],2,glinternet:::standardize,weights=v)
  active <- list(cat=matrix(1,1,1), cont=matrix(1:2,ncol=1),
                 catcat=NULL, contcont=matrix(c(1,2),1,2), catcont=matrix(c(1,1),1,2))
  beta0 <- glinternet:::initial_intercept(y,v,"gaussian")
  beta <- glinternet:::initialize_betahat(active,NULL,beta0,c(2))
  lambda <- fit$lambda[3]
  sol <- glinternet:::group_lasso(Xcat,Z,y,v,active,beta,c(2),lambda,
                                  "gaussian",1e-7,5000,FALSE)
  sizes <- glinternet:::get_group_sizes(sol$activeSet,c(2))
  ends <- cumsum(sizes)
  starts <- c(1,head(ends,-1)+1)
  penalty <- if (length(sizes)) sum(mapply(function(a,b) sqrt(sum(sol$betahat[-1][a:b]^2)),
                                           starts,ends)) else 0
  direct <- sum(v*sol$res^2)/(2*n)+lambda*penalty
  expect_equal(sol$objValue,direct,tolerance=1e-7)
})

test_that("zero-weight contaminated rows do not change fixed-lambda fits", {
  set.seed(72)
  X <- cbind(rnorm(30),rnorm(30))
  y <- 1+X[,1]+rnorm(30,sd=.2)
  Xbad <- rbind(X,c(1e8,-1e8))
  ybad <- c(y,1e12)
  ## n-dependent design scaling is matched by retaining a zero-weight row in both fits.
  Xref <- rbind(X,c(0,0))
  yref <- c(y,0)
  w <- c(rep(1,30),0)
  lambda <- c(.2,.05,.01)
  a <- glinternet(Xbad,ybad,c(1,1),lambda=lambda,weights=w,maxIter=4000)
  b <- glinternet(Xref,yref,c(1,1),lambda=lambda,weights=w,maxIter=4000)
  expect_equal(a$fitted[1:30,],b$fitted[1:30,],tolerance=1e-8)
  expect_equal(a$objValue,b$objValue,tolerance=1e-8)
})

test_that("binomial objective is independent and weight scaling is invariant", {
  set.seed(810)
  n <- 50
  X <- cbind(sample(0:1,n,TRUE),rnorm(n),rnorm(n))
  eta <- -.3+1.1*X[,2]-.7*X[,3]
  y <- rbinom(n,1,plogis(eta))
  w <- runif(n,.2,2)
  v <- glinternet:::validate_weights(w,n)
  Xcat <- matrix(X[,1],ncol=1)
  Z <- apply(X[,2:3,drop=FALSE],2,glinternet:::standardize,weights=v)
  active <- list(cat=matrix(1,1,1),cont=matrix(1:2,ncol=1),catcat=NULL,
                 contcont=matrix(c(1,2),1,2),catcont=matrix(c(1,1),1,2))
  beta <- glinternet:::initialize_betahat(active,NULL,0,c(2))
  sol <- glinternet:::group_lasso(Xcat,Z,y,v,active,beta,c(2),.025,
                                  "binomial",1e-7,6000,FALSE)
  sizes <- glinternet:::get_group_sizes(sol$activeSet,c(2))
  ends <- cumsum(sizes)
  starts <- c(1,head(ends,-1)+1)
  penalty <- if (length(sizes)) sum(mapply(function(a,b) sqrt(sum(sol$betahat[-1][a:b]^2)),
                                           starts,ends)) else 0
  p <- y-sol$res
  etaFit <- qlogis(pmin(pmax(p,1e-15),1-1e-15))
  softplus <- pmax(etaFit,0)+log1p(exp(-abs(etaFit)))
  direct <- sum(v*(softplus-y*etaFit))/n+.025*penalty
  expect_equal(sol$objValue,direct,tolerance=1e-7)

  a <- glinternet(X,y,c(2,1,1),nLambda=4,weights=w,family="binomial",maxIter=4000)
  b <- glinternet(X,y,c(2,1,1),nLambda=4,weights=19*w,family="binomial",maxIter=4000)
  expect_equal(a$lambda,b$lambda,tolerance=1e-12)
  expect_equal(a$fitted,b$fitted,tolerance=1e-5)
})

test_that("final solutions pass full-candidate KKT and represent strong hierarchy", {
  d <- make_data(54)
  fit <- glinternet(d$X,d$y,d$levels,nLambda=6,weights=d$w,
                    maxIter=5000,tol=1e-6)
  i <- length(fit$lambda)
  v <- fit$weights
  Xcat <- d$X[,1:2,drop=FALSE]
  Z <- apply(d$X[,3:4,drop=FALSE],2,glinternet:::standardize,weights=v)
  r <- d$y-fit$fitted[,i]
  candidates <- glinternet:::get_candidates(Xcat,Z,r,v,nrow(d$X),2,2,c(3,2),
                                             NULL,NULL,NULL,NULL)
  checked <- glinternet:::check_kkt(Xcat,Z,r,v,nrow(d$X),2,2,c(3,2),
                                    candidates,fit$activeSet[[i]],fit$lambda[i])
  expect_equal(checked$flag,1)
  expect_equal(checked$activeSet,fit$activeSet[[i]])

  effects <- coef(fit,lambdaIndex=i)[[1]]
  if (!is.null(effects$interactions$catcat))
    expect_true(all(as.vector(effects$interactions$catcat) %in% effects$mainEffects$cat))
  if (!is.null(effects$interactions$contcont))
    expect_true(all(as.vector(effects$interactions$contcont) %in% effects$mainEffects$cont))
  if (!is.null(effects$interactions$catcont)) {
    expect_true(all(effects$interactions$catcont[,1] %in% effects$mainEffects$cat))
    expect_true(all(effects$interactions$catcont[,2] %in% effects$mainEffects$cont))
  }
})

test_that("screened fits satisfy KKT in the retained interaction universe", {
  set.seed(106)
  n <- 60
  X <- matrix(rnorm(n*8),n,8)
  y <- X[,1]+X[,2]+1.2*X[,1]*X[,2]+rnorm(n)
  w <- runif(n)
  fit <- glinternet(X,y,rep(1,8),lambda=.02,weights=w,screenLimit=4,
                    maxIter=5000,tol=1e-6)
  i <- length(fit$lambda)
  v <- fit$weights
  Z <- apply(X,2,glinternet:::standardize,weights=v)
  nullResidual <- y-glinternet:::weighted_mean(y,v)
  retained <- glinternet:::get_candidates(NULL,Z,nullResidual,v,n,0,8,NULL,
                                           NULL,NULL,NULL,4)
  checked <- glinternet:::check_kkt(NULL,Z,y-fit$fitted[,i],v,n,0,8,NULL,
                                    retained,fit$activeSet[[i]],fit$lambda[i])
  expect_true(all(is.finite(fit$objValue)))
  expect_equal(checked$flag,1)
  effects <- coef(fit,lambdaIndex=i)[[1]]
  if (!is.null(effects$interactions$contcont))
    expect_true(all(as.vector(effects$interactions$contcont) %in% effects$mainEffects$cont))
})

test_that("binomial intercept is stationary and frequency replication is equivalent", {
  set.seed(707)
  n <- 36
  X <- cbind(rnorm(n),rnorm(n))
  y <- rbinom(n,1,plogis(-.25+.9*X[,1]-.6*X[,2]))
  w <- sample(1:3,n,replace=TRUE)
  lambda <- .035
  weighted <- glinternet(X,y,c(1,1),lambda=lambda,weights=w,
                          family="binomial",tol=1e-10,maxIter=8000)
  i <- length(weighted$lambda)
  expect_true(abs(sum(weighted$weights*(y-weighted$fitted[,i]))/n) <= 2e-10)

  index <- rep(seq_len(n),w)
  N <- length(index)
  replicated <- glinternet(X[index,,drop=FALSE],y[index],c(1,1),
                            lambda=lambda*sqrt(n/N),family="binomial",
                            tol=1e-10,maxIter=8000)
  expect_equal(predict(weighted,X,"response")[,i],
               predict(replicated,X,"response")[,length(replicated$lambda)],
               tolerance=2e-6)
})
