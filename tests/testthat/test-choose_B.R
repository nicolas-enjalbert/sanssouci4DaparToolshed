test_that("Unit test and fonctional test of choose_B", {
  N1 <- c(2, 100, 10, 5)
  N0 <- c(2, 5, 10, 100)
  for (i in seq_along(N1)) {
    n1 <- N1[i]
    n0 <- N0[i]
    group <- rep(c("1", "0"), times = c(n1, n0))

    expect_true(is.numeric(choose_B(group)))
    expect_equal(choose_B(group), min(1000, choose(n1 + n0, n1)))
  }
})
