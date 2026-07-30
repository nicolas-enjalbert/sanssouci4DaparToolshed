test_that("Unit test for volcanoPlot", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n * p), ncol = n)
  group <- rep(c(1, 0), length.out = n)
  ss_obj <- SansSouci(Y = X, group = group)
  ss_obj <- as.SansSouci4DT(ss_obj)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05
  )
  resplotly <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    interactive = TRUE
  )

  # tester si la sortie est du bon format
  expect_s3_class(res, "ggplot")
  expect_s3_class(resplotly, "plotly")

  # tester si il y a des erreurs en entrée
  ss_obj_fit_nopval <- ss_obj_fit
  ss_obj_fit_nopval$output$p.value <- NULL
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit_nopval,
    pval_thr = 0.5,
    logfc_thr = 0.05
  ))

  ss_obj_fit_noFC <- ss_obj_fit
  ss_obj_fit_noFC$output$estimate <- NULL
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit_noFC,
    pval_thr = 0.5,
    logfc_thr = 0.05
  ))

  ss_obj_fit_noteq <- ss_obj_fit
  ss_obj_fit_noteq$output$estimate <- ss_obj_fit_noteq$output$estimate[-1]
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit_noteq,
    pval_thr = 0.5,
    logfc_thr = 0.05
  ))

  expect_error(res <- volcanoPlot(
    pval_thr = 0.5,
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = matrix(2),
    pval_thr = 0.5,
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = "0.5",
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = c(0.5, 1),
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 1.5,
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = -0.5,
    logfc_thr = 0.05
  ))

  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    qval_thr = "0.5",
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    qval_thr = c(0.5, 1),
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    qval_thr = 1.5,
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    qval_thr = -0.5,
    logfc_thr = 0.05
  ))

  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = "0.05"
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = c(0.05, 1)
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = -0.05
  ))
  expect_warning(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    qval_thr = 0.5,
    logfc_thr = 0.05
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    interactive = "test"
  ))
  expect_error(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    interactive = c(TRUE, TRUE)
  ))
})

test_that("Fonctional test for volcanoPlot", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n * p), ncol = n)
  group <- rep(c(1, 0), length.out = n)
  ss_obj <- SansSouci(Y = X, groups = group)
  ss_obj <- as.SansSouci4DT(ss_obj)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  pval <- as.vector(pValues(ss_obj_fit))
  logFC <- as.vector(foldChanges(ss_obj_fit))
  S <- which((pval < 0.5) & abs(logFC) > 0.05)
  PHB <- predict(ss_obj_fit, S = S)

  res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    condition = c("group1", "group0"),
    pal = c("red", "blue")
  )


  expect_equal(
    res@labels$title,
    paste( #"group1_vs_group0 - ",
      length(S),
      " proteins selected\nAt least ",
      PHB[["TP"]],
      " true positives (FDP \u2264 ",
      sprintf("%.2f", PHB[["FDP"]]),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "Fold change (log scale)")
  expect_equal(res@labels$y, bquote("p-value (-" ~ log[10] ~ "scale)"))

  point_idx <- which(
    sapply(
      res$layers,
      function(x) inherits(x$geom, "GeomPoint")
    )
  )

  point_data <- ggplot_build(res)$data[[point_idx]]

  expect_equal(
    sort(unique(point_data$colour)),
    sort(c("red", "blue"))
  )

  ####### error in condition
  expect_warning(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    conditions = c(
      "group1", "group0",
      "fake_group"
    ),
    pal = c("red", "blue")
  ))


  expect_equal(
    res@labels$title,
    paste( #"group1_vs_group0 - ",
      length(S),
      " proteins selected\nAt least ",
      PHB[["TP"]],
      " true positives (FDP \u2264 ",
      sprintf("%.2f", PHB[["FDP"]]),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "Fold change (log scale)")
  expect_equal(res@labels$y, bquote("p-value (-" ~ log[10] ~ "scale)"))

  point_idx <- which(
    sapply(
      res$layers,
      function(x) inherits(x$geom, "GeomPoint")
    )
  )

  point_data <- ggplot_build(res)$data[[point_idx]]

  expect_equal(
    sort(unique(point_data$colour)),
    sort(c("red", "blue"))
  )

  ####### error in color palette
  expect_warning(res <- volcanoPlot(
    x = ss_obj_fit,
    pval_thr = 0.5,
    logfc_thr = 0.05,
    conditions = c("group1", "group0"),
    pal = c("red", "blue", "grey")
  ))


  expect_equal(
    res@labels$title,
    paste( #"group1_vs_group0 - ",
      length(S),
      " proteins selected\nAt least ",
      PHB[["TP"]],
      " true positives (FDP \u2264 ",
      sprintf("%.2f", PHB[["FDP"]]),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "Fold change (log scale)")
  expect_equal(res@labels$y, bquote("p-value (-" ~ log[10] ~ "scale)"))

  point_idx <- which(
    sapply(
      res$layers,
      function(x) inherits(x$geom, "GeomPoint")
    )
  )

  point_data <- ggplot_build(res)$data[[point_idx]]

  expect_equal(
    sort(unique(point_data$colour)),
    sort(c("orange", "gray"))
  )
})
