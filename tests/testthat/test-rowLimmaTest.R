test_that("unit test of rowLimmaTest() : good use", {
  n <- 5
  p <- 10
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n)

  # for only one design
  res <- rowLimmaTest(Y, groups)

  expect_true(is.list(res))
  expect_equal(length(res), 2)
  expect_equal(names(res), c("p.value", "estimate"))
  expect_equal(length(res$p.value), length(res$estimate))
  expect_equal(length(res$p.value), p)

  # for permuted design
  B <- 7
  null_groups <- replicate(B, sample(groups))
  res <- rowLimmaTest(Y, null_groups)
  expect_true(is.list(res))
  expect_equal(length(res), 2)
  expect_equal(names(res), c("p.value", "estimate"))
  expect_equal(dim(res$p.value), dim(res$estimate))
  expect_equal(dim(res$p.value), c(p, B))
})

test_that("unit test of rowLimmaTest() : bad use", {
  n <- 5
  p <- 10

  ## not same nb of samples
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n-2)

  expect_error(rowLimmaTest(Y, groups), "groups should be of length 5, not 3")

  # for permuted design
  B <- 7
  null_groups <- replicate(B, sample(groups))

  expect_error(rowLimmaTest(Y, null_groups), "groups should be of length 5, not 3")

  ## not good group names
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,2), length.out = n)

  expect_error(rowLimmaTest(Y, groups), regexp =
                 "should consist only of '0' and '1'", fixed = FALSE)
  B <- 7
  null_groups <- replicate(B, sample(groups))
  expect_error(rowLimmaTest(Y, null_groups), regexp =
                 "should consist only of '0' and '1'", fixed = FALSE)

})


test_that("unit test of limmaAnalyses() : good use", {
  n <- 5
  p <- 10
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n)
  res <- limmaAnalyses(Y, groups)

  expect_true(is.list(res))
  expect_equal(length(res), 2)
  expect_equal(names(res), c("p.values", "logFC"))
  expect_equal(length(res$p.values), length(res$logFC))
  expect_equal(length(res$p.values), p)
})

test_that("unit test of limmaAnalyses() : bad use", {
  n <- 5
  p <- 10

  ## not same nb of samples
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n-2)

  expect_error(limmaAnalyses(Y, groups), "groups should be of length 5, not 3")


  ## not good group names
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,2), length.out = n)

  expect_error(limmaAnalyses(Y, groups), regexp =
                 "should consist only of '0' and '1'", fixed = FALSE)

})

test_that("Consistance with Limma analyses", {
  n <- 5
  p <- 10
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n)
  res <- limmaAnalyses(Y, groups)

  res_row <- rowLimmaTest(Y, groups)

  {
    groups <- as.factor(groups)
    # Create design matrix with no intercepts
    design.matrix <- stats::model.matrix(~ 0 + groups)
    # Fit a linear model
    res_lm <- limma::lmFit(Y, design.matrix)
    # Define contrast : here onlye group1 vs group0
    contr <- limma::makeContrasts(groups1 - groups0,
                                  levels = colnames(design.matrix))
    # fit with contrast
    res_fit <- limma::contrasts.fit(res_lm, contr)
    # make test
    res_eb <- limma::eBayes(res_fit)
  }

  expect_equal(as.vector(res$p.values), as.vector(res_eb$p.value))
  expect_equal(as.vector(res$logFC), as.vector(res_eb$coefficients))

  expect_equal(as.vector(res_row$p.value), as.vector(res_eb$p.value))
  expect_equal(as.vector(res_row$estimate), as.vector(res_eb$coefficients))

  B <- 3
  null_groups <- replicate(B, sample(groups))
  res_row_B <- rowLimmaTest(Y, null_groups)
  for(b in seq(B)){
    groupsb <- as.factor(null_groups[,1])
    # Create design matrix with no intercepts
    design.matrix <- stats::model.matrix(~ 0 + groupsb)
    # Fit a linear model
    res_lm <- limma::lmFit(Y, design.matrix)
    # Define contrast : here onlye group1 vs group0
    contr <- limma::makeContrasts(groupsb1 - groupsb0,
                                  levels = colnames(design.matrix))
    # fit with contrast
    res_fit <- limma::contrasts.fit(res_lm, contr)
    # make test
    res_eb <- limma::eBayes(res_fit)

    expect_equal(as.vector(res_row_B$p.value[,1]), as.vector(res_eb$p.value))
    expect_equal(as.vector(res_row_B$estimate[,1]), as.vector(res_eb$coefficients))
  }
})
