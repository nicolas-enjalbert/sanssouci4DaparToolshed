#' @title Extract selection and post hoc bound to QFeature object
#'
#' @description
#' Extract selection and post hoc bound form a `Sanssouci` object to either a
#' `data.frame` a `QFeatures`.
#' This function can be executed multiple times on a given `QFeatures` with
#' different thresholds.
#'
#' @param sanssouci_obj An object of class `SansSouci` which must be calibrated
#' (see fit function in sanssouci package).
#' @param pval_threshold A `numeric(1)`, threshold on pvalues.
#' @param Foldchange_thlogFC A `numeric(1)`, threshold on log Fold Change.
#' @param qfeature_obj An object of class `QFeature`.
#'
#' @returns If qfeature_obj is NULL, returns a `data.frame` containing at least
#' 4 columns :
#' * name_prot : the names of the proteins
#' * pval : the pvalues
#' * logFC : the logFoldChange
#' * sel_n : a binary column for each of the n saved selection
#'
#' Each selection column includes the following attributes:
#' * the number of selected proteins
#' * the thresholds used
#' * the post hoc bounds on the selection.
#'
#' The data.frame includes the following attributes :
#' * the number of permutation used in sanssouci calibration (B)
#' * the test function.
#'
#' If qfeature_obj is not NULL, returns a `QFeature` with the previously
#' described data.frame stored in `metadata(qfeature_obj)$posthoc_selection`.
#'
#' @import sanssouci
#' @import DaparToolshed
#' @importFrom S4Vectors metadata metadata<-
#'
#' @author Nicolas Enjalbert Courrech
#'
#' @export
#'
#' @examples
#' n <- 20
#' p <- 32
#' X <- matrix(rnorm(n * p), ncol = n)
#' group <- rep(c(1, 0), length.out = n)
#' ss_obj <- SansSouci(Y = X, group = group)
#' ss_obj_fit <- fit(ss_obj, alpha = 0.5, B = 100)
#'
#' export2dapartoolshed(ss_obj_fit, 1, 1)
#'
export2dapartoolshed <- function(sanssouci_obj, pval_threshold,
                                 Foldchange_thlogFC, qfeature_obj = NULL) {
  if (missing(sanssouci_obj)) {
    stop("'sanssouci_obj' is required.")
  }
  if (!inherits(sanssouci_obj, "SansSouci")) {
    stop("'sanssouci_obj' must be a 'SansSouci' class object.")
  }

  if (missing(pval_threshold)) {
    stop("'pval_threshold' is required.")
  }
  if (!is.numeric(pval_threshold)) {
    stop("'pval_threshold' must be numeric.")
  }
  if (length(pval_threshold) != 1) {
    stop("'pval_threshold' must be of length 1.")
  }
  if ((pval_threshold < 0) || (pval_threshold > 1)) {
    stop("'pval_threshold' must be between 0 and 1.")
  }
  if (missing(Foldchange_thlogFC)) {
    stop("'Foldchange_thlogFC' is required.")
  }
  if (!is.numeric(Foldchange_thlogFC)) {
    stop("'Foldchange_thlogFC' must be numeric.")
  }
  if (length(Foldchange_thlogFC) != 1) {
    stop("'Foldchange_thlogFC' must be of length 1.")
  }
  if (Foldchange_thlogFC < 0) {
    stop("'Foldchange_thlogFC' must be positive.")
  }

  if (!is.null(qfeature_obj)) {
    if (!inherits(qfeature_obj, "QFeatures")) {
      stop("'qfeature_obj' must be a 'QFeatures' class object.")
    }

    previous_selection <- metadata(qfeature_obj)$posthoc_selection
    if (!is.null(previous_selection)) {
      if (attr(previous_selection, "nb_permutation") != sanssouci_obj$parameters$B) {
        txt <- paste(
          "Number of permutation used must be ",
          attr(previous_selection, "nb_permutation"), "as previously saved."
        )
        stop(txt)
      }
      if (attr(previous_selection, "rowTestFun") != sanssouci_obj$parameters$funName) {
        txt <- paste(
          "Test function used must be ",
          attr(previous_selection, "rowTestFun"), "as previously saved."
        )
        stop(txt)
      }
      df_sel <- previous_selection
    }
  } else {
    previous_selection <- NULL
  }

  pval <- sanssouci::pValues(sanssouci_obj)
  logFC <- sanssouci::foldChanges(sanssouci_obj)


  sel_prot <- (pval < pval_threshold) & (abs(logFC) > Foldchange_thlogFC)

  names_prot <- rownames(sanssouci_obj$input$Y)
  if (is.null(names_prot)) {
    names_prot <- seq_along(pval)
  }
  # name_sel_prot <- names_prot[sel_prot]

  if (is.null(previous_selection)) {
    df_sel <- data.frame(
      name_prot = names_prot,
      pval = as.vector(pval),
      logFC = as.vector(logFC)
    )
    attr(df_sel, "nb_permutation") <- sanssouci_obj$parameters$B
    attr(df_sel, "rowTestFun") <- sanssouci_obj$parameters$funName
  }
  pred_sel <- predict(sanssouci_obj, S = which(sel_prot))
  name_col_sel <- paste("sel", ncol(df_sel) - 2, sep = "_")
  df_sel[[name_col_sel]] <- as.vector(sel_prot * 1)
  attributes(df_sel[[name_col_sel]]) <- list(
    thr_pval = pval_threshold,
    thr_logFC = Foldchange_thlogFC,
    n_sel = sum(sel_prot),
    FDP = pred_sel[["FDP"]],
    TP = pred_sel[["TP"]]
  )

  if (is.null(qfeature_obj)) {
    return(df_sel)
  } else {
    metadata(qfeature_obj)$posthoc_selection <- df_sel

    return(qfeature_obj)
  }
}
