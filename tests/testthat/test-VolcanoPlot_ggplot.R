test_that("Unit test for VolcanoPlot_ss4DT", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n * p), ncol = n)
  group <- rep(c(1, 0), length.out = n)
  ss_obj <- SansSouci(Y = X, group = group)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05
  )
  resplotly <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    interactive = TRUE
  )

  # tester si la sortie est du bon format
  expect_s3_class(res, "ggplot")
  expect_s3_class(resplotly, "plotly")

  # tester si il y a des erreurs en entrée
  ss_obj_fit_nopval <- ss_obj_fit
  ss_obj_fit_nopval$output$p.value <- NULL
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit_nopval,
    th_pval = 0.5,
    th_logfc = 0.05
  ))

  ss_obj_fit_noFC <- ss_obj_fit
  ss_obj_fit_noFC$output$estimate <- NULL
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit_noFC,
    th_pval = 0.5,
    th_logfc = 0.05
  ))

  ss_obj_fit_noteq <- ss_obj_fit
  ss_obj_fit_noteq$output$estimate <- ss_obj_fit_noteq$output$estimate[-1]
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit_noteq,
    th_pval = 0.5,
    th_logfc = 0.05
  ))

  expect_error(res <- VolcanoPlot_ss4DT(
    th_pval = 0.5,
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = matrix(2),
    th_pval = 0.5,
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = "0.5",
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = c(0.5, 1),
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 1.5,
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = -0.5,
    th_logfc = 0.05
  ))

  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_qval = "0.5",
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_qval = c(0.5, 1),
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_qval = 1.5,
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_qval = -0.5,
    th_logfc = 0.05
  ))

  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = "0.05"
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = c(0.05, 1)
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = -0.05
  ))
  expect_warning(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_qval = 0.5,
    th_logfc = 0.05
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    interactive = "test"
  ))
  expect_error(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    interactive = c(TRUE, TRUE)
  ))
})

test_that("Fonctional test for VolcanoPlot_ss4DT", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n * p), ncol = n)
  group <- rep(c(1, 0), length.out = n)
  ss_obj <- SansSouci(Y = X, group = group)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  pval <- as.vector(pValues(ss_obj_fit))
  logFC <- as.vector(foldChanges(ss_obj_fit))
  S <- which((pval < 0.5) & abs(logFC) > 0.05)
  PHB <- predict(ss_obj_fit, S = S)

  res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    condition = c("group1", "group0"),
    pal = c("red", "blue")
  )


  expect_equal(
    res@labels$title,
    paste("group1_vs_group0 - ",
      length(S),
      " proteins selected \nAt least ",
      PHB[["TP"]],
      " true positives (FDP <= ",
      round(PHB[["FDP"]], 2),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "logFC")
  expect_equal(res@labels$y, "-log10(pValue)")

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
  expect_warning(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    conditions = c(
      "group1", "group0",
      "fake_group"
    ),
    pal = c("red", "blue")
  ))


  expect_equal(
    res@labels$title,
    paste(length(S),
      " proteins selected \nAt least ",
      PHB[["TP"]],
      " true positives (FDP <= ",
      round(PHB[["FDP"]], 2),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "logFC")
  expect_equal(res@labels$y, "-log10(pValue)")

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
  expect_warning(res <- VolcanoPlot_ss4DT(
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05,
    conditions = c("group1", "group0"),
    pal = c("red", "blue", "grey")
  ))


  expect_equal(
    res@labels$title,
    paste("group1_vs_group0 - ",
      length(S),
      " proteins selected \nAt least ",
      PHB[["TP"]],
      " true positives (FDP <= ",
      round(PHB[["FDP"]], 2),
      ")",
      sep = ""
    )
  )
  expect_equal(res@labels$x, "logFC")
  expect_equal(res@labels$y, "-log10(pValue)")

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
