test_that("Unit test for VolcanoPlot_ggplot", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n*p), ncol = n)
  group <- rep(c(1,0), length.out = n)
  ss_obj <- SansSouci(Y = X, group = group)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  pval <- as.vector(pValues(ss_obj_fit))
  logFC <- as.vector(foldChanges(ss_obj_fit))

  res <- VolcanoPlot_ggplot(pval = pval,
                            logFC = logFC,
                            sanssouci_object = ss_obj_fit,
                            th_pval = 0.5,
                            th_logfc = 0.05)

  # tester si la sortie est du bon format
  expect_s3_class(res, "ggplot")

  # tester si il y a des erreurs en entrée
  expect_error(res <- VolcanoPlot_ggplot(
    logFC = logFC,
    sanssouci_object = ss_obj_fit,
    th_pval = 0.5,
    th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = as.character(pval),
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,

                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = as.character(logFC),
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = 1,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,

                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = matrix(2),
                                         th_pval = 0.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = "0.5",
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = c(0.5, 1),
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 1.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = -0.5,
                                         th_logfc = 0.05))

  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_qval = "0.5",
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_qval = c(0.5, 1),
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_qval = 1.5,
                                         th_logfc = 0.05))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_qval = -0.5,
                                         th_logfc = 0.05))

  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = "0.05"))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = c(0.05, 1)))
  expect_error(res <- VolcanoPlot_ggplot(pval = pval,
                                         logFC = logFC,
                                         sanssouci_object = ss_obj_fit,
                                         th_pval = 0.5,
                                         th_logfc = -0.05))
  expect_warning(res <- VolcanoPlot_ggplot(pval = pval,
                                           logFC = logFC,
                                           sanssouci_object = ss_obj_fit,
                                           th_pval = 0.5,
                                           th_qval = 0.5,
                                           th_logfc = 0.05))


})

test_that("Fonctional test for VolcanoPlot_ggplot", {
  n <- 20
  p <- 32
  X <- matrix(rnorm(n*p), ncol = n)
  group <- rep(c(1,0), length.out = n)
  ss_obj <- SansSouci(Y = X, group = group)
  ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)

  pval <- as.vector(pValues(ss_obj_fit))
  logFC <- as.vector(foldChanges(ss_obj_fit))
  S = which((pval < 0.5) & abs(logFC) > 0.05)
  PHB <- predict(ss_obj_fit, S = S)

  res <- VolcanoPlot_ggplot(pval = pval,
                            logFC = logFC,
                            sanssouci_object = ss_obj_fit,
                            th_pval = 0.5,
                            th_logfc = 0.05,
                            condition = c("group1", "group0"),
                            pal = c("red", "blue"))


  expect_equal(res@labels$title,
               paste("group1_vs_group0 - ",
                     length(S),
                     " proteins selected \nAt least ",
                     PHB[["TP"]],
                     " true positives (FDP <= ",
                     round(PHB[["FDP"]], 2),
                     ")", sep = ""))
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
  expect_warning(res <- VolcanoPlot_ggplot(pval = pval,
                                           logFC = logFC,
                                           sanssouci_object = ss_obj_fit,
                                           th_pval = 0.5,
                                           th_logfc = 0.05,
                                           conditions = c("group1", "group0",
                                                          "fake_group"),
                                           pal = c("red", "blue")))


  expect_equal(res@labels$title,
               paste(length(S),
                     " proteins selected \nAt least ",
                     PHB[["TP"]],
                     " true positives (FDP <= ",
                     round(PHB[["FDP"]], 2),
                     ")", sep = ""))
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
  expect_warning(res <- VolcanoPlot_ggplot(pval = pval,
                                           logFC = logFC,
                                           sanssouci_object = ss_obj_fit,
                                           th_pval = 0.5,
                                           th_logfc = 0.05,
                                           conditions = c("group1", "group0"),
                                           pal = c("red", "blue", "grey")))


  expect_equal(res@labels$title,
               paste("group1_vs_group0 - ",
                     length(S),
                     " proteins selected \nAt least ",
                     PHB[["TP"]],
                     " true positives (FDP <= ",
                     round(PHB[["FDP"]], 2),
                     ")", sep = ""))
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
