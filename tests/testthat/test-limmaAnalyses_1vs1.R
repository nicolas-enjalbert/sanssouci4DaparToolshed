test_that("unit test of limmaAnalyses_1vs1()", {
  n <- 5
  p <- 10
  Y <- matrix(rnorm(n*p), ncol = n)
  groups = rep_len(c(1,0), length.out = n)
  res <- limmaAnalyses_1vs1(Y, groups)

  expect_true(is.list(res))
  expect_equal(length(res), 2)
  expect_equal(names(res), c("p.values", "logFC"))
  expect_equal(length(res$p.values), length(res$logFC))
  expect_equal(length(res$p.values), p)
})
