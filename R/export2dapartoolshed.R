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
#' @param pval_thr A `numeric(1)`, threshold on pvalues.
#' @param logfc_thr A `numeric(1)`, threshold on log Fold Change.
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
#' @importFrom stats predict
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
#' ss_obj <- sanssouci::SansSouci(Y = X, group = group)
#' ss_obj_fit <- sanssouci::fit(ss_obj, alpha = 0.5, B = 100)
#'
#' export2dapartoolshed(ss_obj_fit, 1, 1)
#'
export2dapartoolshed <- function(sanssouci_obj, pval_thr,
                                 logfc_thr, qfeature_obj = NULL) {
  if (missing(sanssouci_obj)) {
    stop("'sanssouci_obj' is required.")
  }
  if (!inherits(sanssouci_obj, "SansSouci")) {
    stop("'sanssouci_obj' must be a 'SansSouci' class object.")
  }

  if (missing(pval_thr)) {
    stop("'pval_thr' is required.")
  }
  if (!is.numeric(pval_thr)) {
    stop("'pval_thr' must be numeric.")
  }
  if (length(pval_thr) != 1) {
    stop("'pval_thr' must be of length 1.")
  }
  if ((pval_thr < 0) || (pval_thr > 1)) {
    stop("'pval_thr' must be between 0 and 1.")
  }
  if (missing(logfc_thr)) {
    stop("'logfc_thr' is required.")
  }
  if (!is.numeric(logfc_thr)) {
    stop("'logfc_thr' must be numeric.")
  }
  if (length(logfc_thr) != 1) {
    stop("'logfc_thr' must be of length 1.")
  }
  if (logfc_thr < 0) {
    stop("'logfc_thr' must be positive.")
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


  sel_prot <- (pval < pval_thr) & (abs(logFC) > logfc_thr)

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
    pval_thr = pval_thr,
    logfc_thr = logfc_thr,
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




#' Get info of post hoc selection
#'
#' @param object QFeatures object containing post hoc selection
#' @param selection_name character, name of the selection
#'
#' @importFrom S4Vectors metadata metadata<-
#'
#' @returns a list containing:
#' * pval_thr: the threshold on pvalues used for the selection,
#' * logfc_thr: the threshold on log Foldchange used for the selection,
#' * n_sel: the number of selected proteins,
#' * FDP: the post hoc bound on the FDP,
#' * TP: the post hoc bound on the TP
#'
#' @author Nicolas Enjalbert-Courrech and Pierre Neuvial
#'
#' @export
#'
#' @examples
#' numdata <- data.frame(
#'   S1 = c(2, 5, 9, 10, 7, 8, 6, 8, 6, 7),
#'   S2 = c(7, 8, 5, 5, 7, 9, 9, 4, 9, 1),
#'   S3 = c(4, 1, 4, 9, 6, 8, 5, 3, 4, 3),
#'   S4 = c(9, 8, 9, 8, 3, 6, 4, 1, 4, 5),
#'   names = paste0("Prot", seq_len(10))
#' )
#' obj <- QFeatures::readQFeatures(numdata,
#'   quantCols = seq_len(4),
#'   fnames = "names",
#'   name = "datatest"
#' )
#'
#' coldata <- data.frame(
#'   quantCols = c("S1", "S2", "S3", "S4"),
#'   Condition = c("C1", "C1", "C2", "C2"),
#'   Bio.Rep = as.character(seq_len(4))
#' )
#' rownames(coldata) <- coldata$quantCols
#' SummarizedExperiment::colData(obj) <- coldata
#'
#' ss_obj <- SansSouci4DT(obj)
#' ss_obj <- sanssouci::fit(ss_obj, alpha = 0.5, B = 0)
#' obj <- export2dapartoolshed(ss_obj, 0.1, 0.5, obj)
#' getPostHocBound(obj, "sel_1")
getPostHocBound <- function(object, selection_name){

  if (missing(object)) {
    stop("'object' is required.")
  }
  if (missing(selection_name)) {
    stop("'selection_name' is required.")
  }
  if(!inherits(object, "QFeatures")){
    stop("object must be a QFeatures object.")
  }

  df <- metadata(object)

  if(is.null(df$posthoc_selection)){
    stop("object must contain 'posthoc_selection' data.frame.")
  }
  if(is.null(df$posthoc_selection[[selection_name]])){
    stop(selection_name, " is not save in posthoc_selection data.frame in object.")
  }

  attributes(df$posthoc_selection[[selection_name]])
}
