#' @title Convert a QFeatures in Sanssouci object
#'
#' @description 
#' Convert an assay of an object of class `QFeatures` in a `Sanssouci` object.
#'
#' @param obj An object of class `QFeatures`.
#' @param i An `integer(1)` index or a `character(1)` name of the assay which will be converted. 
#' If NULL, the last assay will be selected by default.
#'
#' @returns An object of class `SansSouci`.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom DaparToolshed design_qf
#' @importFrom sanssouci SansSouci
#'
#' @author Manon Gaudin
#'
#' @export
#'
#' @examples
#' #xxx a faire
#' NULL
#' 

wrapper_QFtoSansSouci <- function(obj, i = NULL) {
  if (missing(obj)){
    stop("'obj' is required.")
  }
  if (!is(obj, "QFeatures")){
    stop("'obj' must be an object of class QFeatures.")
  }
  if (is.null(i)){
    i <- length(obj)
  }
  if (length(i) != 1){
    stop("'i' must be of length 1.")
  }
  if (!is.numeric(i) & !(i %in% names(obj))){
    stop("'i' must either be numeric or be the name of one obj assay.")
  }
  if (is.numeric(i) & (i > length(obj))){
    txt <- paste0("'i' is out of bounds. Numeric 'i' must be between 1 and ",
                  length(obj), ".")
    stop(txt)
  }
  
  Y <- SummarizedExperiment::assay(obj[[i]])
  groups <- as.integer(as.factor(DaparToolshed::design_qf(obj)$Condition)) - 1 
  
  ## gerer le cas où + que 2 conditions ?
  
  SansSouciobj <- sanssouci::SansSouci(Y = Y, 
                                       groups = groups)
  return(SansSouciobj)
}
