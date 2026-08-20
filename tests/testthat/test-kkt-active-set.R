test_that("KKT updates retain every added one-column main-effect row", {
  Z <- cbind(
    c(1, 0, 0, -1),
    c(0, 1, 0, -1),
    c(0, 0, 1, -1)
  )
  residual <- c(1, 1, 1, -3)
  group_names <- c("cat", "cont", "catcat", "contcont", "catcont")

  variables <- setNames(vector("list", 5), group_names)
  variables$cont <- matrix(1:3, ncol=1)
  candidates <- list(
    variables=variables,
    norms=setNames(vector("list", 5), group_names)
  )
  active_set <- setNames(vector("list", 5), group_names)
  active_set$cont <- matrix(1, ncol=1)

  expect_no_warning(
    checked <- glinternet:::check_kkt(
      X=NULL, Z=Z, res=residual, weights=rep(1, 4), n=4,
      pCat=0, pCont=3, numLevels=NULL, candidates=candidates,
      activeSet=active_set, lambda=.1
    )
  )

  expect_equal(checked$flag, 0)
  expect_equal(checked$activeSet$cont, matrix(1:3, ncol=1))
})
