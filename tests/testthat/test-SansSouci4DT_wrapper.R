test_that("unit test of SansSouci4DT()", {
  numdata <- data.frame(
    S1 = c(2, 5, 9, 10, 7, 8, 6, 8, 6, 7),
    S2 = c(7, 8, 5, 5, 7, 9, 9, 4, 9, 1),
    S3 = c(4, 1, 4, 9, 6, 8, 5, 3, 4, 3),
    S4 = c(9, 8, 9, 8, 3, 6, 4, 1, 4, 5),
    names = paste0("Prot", seq_len(10))
  )
  obj <- QFeatures::readQFeatures(numdata, quantCols = seq_len(4), fnames = "names", name = "datatest")

  coldata <- data.frame(
    quantCols = c("S1", "S2", "S3", "S4"),
    Condition = c("C1", "C1", "C2", "C2"),
    Bio.Rep = as.character(seq_len(4))
  )
  rownames(coldata) <- coldata$quantCols
  SummarizedExperiment::colData(obj) <- coldata

  res <- SansSouci4DT(obj)

  expect_is(res, "SansSouci")
  expect_is(res, "SansSouci4DT")
  expect_equal(length(res), 3)
  expect_true(all(res$input$Y == numdata[, -5]))

  expect_error(SansSouci4DT())
  expect_error(SansSouci4DT(obj[[1]]))
  expect_error(SansSouci4DT(obj, i = c(1, 2)))
  expect_error(SansSouci4DT(obj, i = "test"))
  expect_error(SansSouci4DT(obj, i = 3))
})
