#' @title Convert a QFeatures in Sanssouci object
#'
#' @description
#' Convert an assay of an object of class `QFeatures` in a `Sanssouci` object.
#'
#' @param obj An object of class `QFeatures`.
#' @param i An `integer(1)` index or a `character(1)` name of the assay which
#' will be converted. If NULL, the last assay will be selected by default.
#'
#' @returns An object of class `SansSouci4DT`.
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
#' wrapper_QFtoSansSouci(obj)
#'
wrapper_QFtoSansSouci <- function(obj, i = NULL) {
  if (missing(obj)) {
    stop("'obj' is required.")
  }
  if (!inherits(obj, "QFeatures")) {
    stop("'obj' must be an object of class QFeatures.")
  }
  if (is.null(i)) {
    i <- length(obj)
  }
  if (length(i) != 1) {
    stop("'i' must be of length 1.")
  }
  if (!is.numeric(i) && !(i %in% names(obj))) {
    stop("'i' must either be numeric or be the name of one obj assay.")
  }
  if (is.numeric(i) && (i > length(obj))) {
    txt <- paste0(
      "'i' is out of bounds. Numeric 'i' must be between 1 and ",
      length(obj), "."
    )
    stop(txt)
  }

  Y <- SummarizedExperiment::assay(obj[[i]])
  groups <- as.integer(as.factor(DaparToolshed::design_qf(obj)$Condition)) - 1

  ## gerer le cas où + que 2 conditions ?

  SansSouciobj <- sanssouci::SansSouci(
    Y = Y,
    groups = groups
  )
  SansSouciobj <- as.SansSouci4DT(SansSouciobj)
  return(SansSouciobj)
}
