test_that("Sanity checks of 'categCheck' throw errors when expected to", {
  n <- 10

  vec <- rep(1, 5)
  expect_error(
    categCheck(vec, n),
    "vec should be of length 10, not 5"
  )

  expect_no_error(categCheck(c(0, 1, 0, 1), 4))

  expect_error(categCheck(c(1, 2, 1, 2), 4),
    regexp =
      "should consist only of '0' and '1'", fixed = FALSE
  )

  expect_no_error(categCheck(c(1, 2, 3, 4), 4))
})

test_that("UT for as.SansSouci4DT", {
  ss_obj <- sanssouci::SansSouciSim(m = 10, pi0 = 1, n = 5, rho = 1)
  ss4DT_obj <- as.SansSouci4DT(ss_obj)
  expect_equal(class(ss4DT_obj), c("SansSouci4DT", "SansSouci"))
  expect_error(as.SansSouci4DT(1))
})
