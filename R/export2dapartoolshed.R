#' Extract selection and post hoc bound to QFeature object
#'
#' @param sanssouci_obj object of class 'SansSouci' which must be calibrated
#' (see fit function in sanssouci package)
#' @param pval_threshold numeric, threshold on pvalues
#' @param Foldchange_thlogFC numeric, threshold on log Fold Change
#' @param qfeature_obj object of class 'QFeature'. If NULL, this function
#' returns a data.frame. If not, the data.frame is stored in
#' metadata(qfeature_obj)$posthoc_selection
#'
#' @returns a data.frame with the names of the proteins (name_prot),
#' the pvalues (pval), the logfoldchange (logFC) and binary columns for each
#' saved selection. Each selection column includes the following attributes:
#' the number of selected proteins, the thresholds used, and the post hoc bounds
#' on the selection. The data.frame includes the following attributes :
#' the number of permutation used in sanssouci calibration (B) and the test
#' function.
#' If qfeature_objec os NULL, this function
#' returns a data.frame. If not, the data.frame is stored in
#' metadata(qfeature_obj)$posthoc_selection
#'
#' @export
#'
#' @examples
#' #xxx a faire
#' NULL
export2dapartoolshed <- function(sanssouci_obj, pval_threshold,
                                 Foldchange_thlogFC, qfeature_obj = NULL){
  if(class(sanssouci_obj) != "SansSouci") {
    stop("sanssouci_obj must be a 'SansSouci' class object.")
  }

  if(!is.null(qfeature_obj)){
  if(class(qfeature_obj) != "QFeatures") {
    stop("qfeature_obj must be a 'QFeatures' class object.")
  }

  previous_selection <- metadata(qfeature_obj)$posthoc_selection
  if(!is.null(previous_selection)){
    if(attr(previous_selection, "nb_permutation") != sanssouci_obj$parameters$B){
      stop(paste("Number of permutation used must be ",
                 attr(previous_selection, "nb_permutation"), "as previously saved.")
      )
    }
    if(attr(previous_selection, "rowTestFun") != sanssouci_obj$parameters$funName){
      stop(paste("Test function used must be ",
                 attr(previous_selection, "rowTestFun"), "as previously saved.")
      )
    }
    df_sel <- previous_selection
  }

  } else {
    previous_selection <- NULL
  }

  pval <- pValues(sanssouci_obj)
  logFC <- foldChanges(sanssouci_obj)


  sel_prot <- (pval < pval_threshold) & (abs(logFC) > Foldchange_thlogFC)

  names_prot <- rownames(sanssouci_obj$input$Y)

  name_sel_prot <- names_prot[sel_prot]

  if(is.null(previous_selection)){
    df_sel <- data.frame(name_prot = names_prot,
                         pval = as.vector(pval),
                         logFC = as.vector(logFC))
    attr(df_sel, "nb_permutation") <- sanssouci_obj$parameters$B
    attr(df_sel, "rowTestFun") <- sanssouci_obj$parameters$funName
  }
  pred_sel <- predict(sanssouci_obj, S = which(sel_prot))
  name_col_sel <- paste("sel", ncol(df_sel)-2, sep = "_")
  df_sel[[name_col_sel]] <- as.vector(sel_prot*1)
  attributes(df_sel[[name_col_sel]]) <- list(thr_pval = pval_threshold,
                                             thr_logFC = Foldchange_thlogFC,
                                             n_sel = sum(sel_prot),
                                             FDP = pred_sel[["FDP"]],
                                             TP = pred_sel[["TP"]])

  if(is.null(qfeature_obj)){
    return(df_sel)
  } else {
  metadata(qfeature_obj)$posthoc_selection <- df_sel

  return(qfeature_obj)
  }
}



