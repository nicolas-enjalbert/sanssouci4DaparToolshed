library(sanssouci)

numdata <- data.frame(
  S1 = c(2, 5, 9, 10, 7, 8, 6, 8, 6, 7),
  S2 = c(7, 8, 5, 5, 7, 9, 9, 4, 9, 1),
  S3 = c(4, 1, 4, 9, 6, 8, 5, 3, 4, 3),
  S4 = c(9, 8, 9, 8, 3, 6, 4, 1, 4, 5),
  names = paste0("Prot", seq_len(10))
)
nb_prot <- nrow(numdata)
nb_samples <- ncol(numdata) - 1
qfeature_obj <- QFeatures::readQFeatures(numdata, quantCols = seq_len(4),
                                         fnames = "names", name = "datatest")

coldata <- data.frame(
  quantCols = c("S1", "S2", "S3", "S4"),
  Condition = c("C1", "C1", "C2", "C2"),
  Bio.Rep = as.character(seq_len(4))
)
rownames(coldata) <- coldata$quantCols
SummarizedExperiment::colData(qfeature_obj) <- coldata

# rownames(numdata) <- numdata$names
# numdata <- numdata[, colnames(numdata) != "names"]
sanssouci_init <- SansSouci4DT(qfeature_obj)

B <- 0
sanssouci_obj <- fit(sanssouci_init, B = B, alpha = 0.5)


test_that("Unit/fctal test of export2dapartoolshed : good", {
  ## dont give qfeature_obj ; expect data.frame object
  pval_thr <- 0.5
  logFC_thr <- 0.001
  res <- export2dapartoolshed(sanssouci_obj,
    pval_thr = pval_thr,
    logfc_thr = logFC_thr
  )

  expect_true(is.data.frame(res))
  expect_equal(colnames(res), c("name_prot", "pval", "logFC", "sel_1"))
  expect_equal(names(attributes(res)), c(
    "names", "row.names", "nb_permutation",
    "rowTestFun", "class"
  ))

  sanssouci_obj_tmp <- sanssouci_obj
  rownames(sanssouci_obj_tmp$input$Y) <- NULL
  res2 <- export2dapartoolshed(sanssouci_obj_tmp,
                              pval_thr = pval_thr,
                              logfc_thr = logFC_thr
  )

  expect_equal(res2$name_prot, seq_len(10))

  # expect_equal(names(attributes(res$sel_1)), c("thr_pval", "thr_logFC",
  #                                              "n_sel",
  #                                       "FDP", "TP")) # tested with consistenc
  expect_equal(dim(res), c(nb_prot, 4))
  expect_equal(res$name_prot, numdata$names)

  ## Consistency with sanssouci
  pval <- pValues(sanssouci_obj)
  logFC <- foldChanges(sanssouci_obj)
  S <- (pval < pval_thr) & (abs(logFC) > logFC_thr)
  PHB <- predict(sanssouci_obj, S = which(S), what = c("FDP", "TP"))

  expect_equal(as.vector(res$pval), as.vector(pval))
  expect_equal(as.vector(res$logFC), as.vector(logFC))
  expect_equal(as.vector(res$sel_1), as.vector(S * 1))

  expect_equal(attr(res, "nb_permutation"), B)
  expect_equal(attr(res, "rowTestFun"), "rowWelchTests")

  expect_equal(
    attributes(res$sel_1),
    list(
      "pval_thr" = pval_thr,
      "logfc_thr" = logFC_thr,
      "n_sel" = sum(S),
      "FDP" = PHB[["FDP"]],
      "TP" = PHB[["TP"]]
    )
  )

  #### test with qfeature_obj : expect QFeature object


  res1 <- export2dapartoolshed(
    qfeature_obj = qfeature_obj,
    sanssouci_obj = sanssouci_obj,
    pval_thr = pval_thr,
    logfc_thr = logFC_thr
  )
  expect_true(inherits(res1, class(qfeature_obj)))
  res <- metadata(res1)$posthoc_selection
  expect_true(is.data.frame(res))
  expect_equal(colnames(res), c("name_prot", "pval", "logFC", "sel_1"))
  expect_equal(names(attributes(res)), c(
    "names", "row.names", "nb_permutation",
    "rowTestFun", "class"
  ))
  # expect_equal(names(attributes(res$sel_1)), c("thr_pval", "thr_logFC",
  #                                              "n_sel",
  #                                       "FDP", "TP")) # tested with consistenc
  expect_equal(dim(res), c(nb_prot, 4))
  expect_equal(res$name_prot, numdata$names)

  ## Consistency with sanssouci
  pval <- pValues(sanssouci_obj)
  logFC <- foldChanges(sanssouci_obj)
  S <- (pval < pval_thr) & (abs(logFC) > logFC_thr)
  PHB <- predict(sanssouci_obj, S = which(S), what = c("FDP", "TP"))

  expect_equal(as.vector(res$pval), as.vector(pval))
  expect_equal(as.vector(res$logFC), as.vector(logFC))
  expect_equal(as.vector(res$sel_1), as.vector(S * 1))

  expect_equal(attr(res, "nb_permutation"), B)
  expect_equal(attr(res, "rowTestFun"), "rowWelchTests")

  expect_equal(
    attributes(res$sel_1),
    list(
      "pval_thr" = pval_thr,
      "logfc_thr" = logFC_thr,
      "n_sel" = sum(S),
      "FDP" = PHB[["FDP"]],
      "TP" = PHB[["TP"]]
    )
  )

  pval_thr2 <- 0.05
  logFC_thr2 <- 0.04
  res2 <- export2dapartoolshed(
    qfeature_obj = res1,
    sanssouci_obj = sanssouci_obj,
    pval_thr = pval_thr2,
    logfc_thr = logFC_thr2
  )
  expect_true(inherits(res2, class(qfeature_obj)))
  res <- metadata(res2)$posthoc_selection
  expect_true(is.data.frame(res))
  expect_equal(colnames(res), c("name_prot", "pval", "logFC", "sel_1", "sel_2"))
  expect_equal(names(attributes(res)), c(
    "names", "row.names", "nb_permutation",
    "rowTestFun", "class"
  ))
  # expect_equal(names(attributes(res$sel_1)), c("thr_pval", "thr_logFC",
  #                                              "n_sel",
  #                                       "FDP", "TP")) # tested with consistenc
  expect_equal(dim(res), c(nb_prot, 5))
  expect_equal(res$name_prot, numdata$names)

  ## Consistency with sanssouci : sel_1 is not modified
  pval <- pValues(sanssouci_obj)
  logFC <- foldChanges(sanssouci_obj)
  S <- (pval < pval_thr) & (abs(logFC) > logFC_thr)
  PHB <- predict(sanssouci_obj, S = which(S), what = c("FDP", "TP"))

  expect_equal(as.vector(res$pval), as.vector(pval))
  expect_equal(as.vector(res$logFC), as.vector(logFC))
  expect_equal(as.vector(res$sel_1), as.vector(S * 1))

  expect_equal(attr(res, "nb_permutation"), B)
  expect_equal(attr(res, "rowTestFun"), "rowWelchTests")

  expect_equal(
    attributes(res$sel_1),
    list(
      "pval_thr" = pval_thr,
      "logfc_thr" = logFC_thr,
      "n_sel" = sum(S),
      "FDP" = PHB[["FDP"]],
      "TP" = PHB[["TP"]]
    )
  )

  ## Consistency with sanssouci : test sel_2
  S <- (pval < pval_thr2) & (abs(logFC) > logFC_thr2)
  PHB <- predict(sanssouci_obj, S = which(S), what = c("FDP", "TP"))

  expect_equal(as.vector(res$sel_2), as.vector(S * 1))

  expect_equal(
    attributes(res$sel_2),
    list(
      "pval_thr" = pval_thr2,
      "logfc_thr" = logFC_thr2,
      "n_sel" = sum(S),
      "FDP" = PHB[["FDP"]],
      "TP" = PHB[["TP"]]
    )
  )
})

test_that("Unit test of export2dapartoolshed : errors", {
  # Errors with sanssouci_obj argument
  expect_error(export2dapartoolshed(
      pval_thr = 1,
      logfc_thr = 1
    ))
  # Do not give a sanssouci object
  expect_error(
    export2dapartoolshed(
      sanssouci_obj = 10,
      pval_thr = 1,
      logfc_thr = 1
    )
  )

  # Errors with pval_thr argument
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    logfc_thr = 1
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = "1",
    logfc_thr = 1
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = c(1,1),
    logfc_thr = 1
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = -1,
    logfc_thr = 1
  ))

  # Errors with logfc_thr argument
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = 1,
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = 1,
    logfc_thr = "1"
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = 1,
    logfc_thr = c(1,1)
  ))
  expect_error(export2dapartoolshed(
    sanssouci_obj = sanssouci_obj,
    pval_thr = 1,
    logfc_thr = -1
  ))

  # Do not give a qfeatures object
  expect_error(
    export2dapartoolshed(
      sanssouci_obj = sanssouci_obj,
      qfeature_obj = 10,
      pval_thr = 1,
      logfc_thr = 1
    ),
    regexp = "'qfeature_obj' must be a 'QFeatures' class object."
  )


  res1 <- export2dapartoolshed(
    qfeature_obj = qfeature_obj,
    sanssouci_obj = sanssouci_obj,
    pval_thr = 0.5,
    logfc_thr = 0.5
  )
  # try to use another calibration parameters : B
  sanssouci_obj_F <- fit(sanssouci_init, alpha = 0.5, B = 1)
  expect_error(
    export2dapartoolshed(
      sanssouci_obj = sanssouci_obj_F,
      pval_thr = 1,
      qfeature_obj = res1,
      logfc_thr = 1
    ),
    regexp = "Number of permutation used must be  0 as previously saved."
  )
  # try to use another calibration parameters : rowTestFun
  sanssouci_obj_F <- fit(sanssouci_init,
    alpha = 0.5,
    rowTestFUN = rowWilcoxonTests, B = 0
  )
  expect_error(
    export2dapartoolshed(
      sanssouci_obj = sanssouci_obj_F,
      pval_thr = 1,
      qfeature_obj = res1,
      logfc_thr = 1
    ),
    regexp = "Test function used must be  rowWelchTests as previously saved."
  )
})

test_that("Unit test of getPostHocBound", {
  pval_thr <- 0.5
  logFC_thr <- 0.001
  qf_obj_PH <- export2dapartoolshed(sanssouci_obj,
                              pval_thr = pval_thr,
                              logfc_thr = logFC_thr,
                              qfeature_obj = qfeature_obj
  )

  res <- getPostHocBound(qf_obj_PH, selection_name = "sel_1")
  expect_is(res, "list")
  expect_equal(names(res), c("pval_thr", "logfc_thr", "n_sel", "FDP", "TP"))
  expect_equal(res$pval_thr, pval_thr)
  expect_equal(res$logfc_thr, logFC_thr)
  expect_equal(res$n_sel, 5)
  expect_equal(res$FDP, 1)
  expect_equal(res$TP, 0)

  ## error
  expect_error({getPostHocBound(selection_name = "sel_1")},
              "'object' is required")
  expect_error({getPostHocBound(object = matrix(NA, nrow = 3, ncol = 4),
                                selection = "sel_1")},
               "object must be a QFeatures object.")
  expect_error({getPostHocBound(object = qf_obj_PH)},
               "'selection_name' is required")
  expect_error({getPostHocBound(object = qfeature_obj, selection = "sel_1")},
               "object must contain 'posthoc_selection' data.frame.")
  expect_error({getPostHocBound(object = qf_obj_PH, selection = "wrongname")},
               "wrongname is not save in posthoc_selection data.frame in object.")
})
