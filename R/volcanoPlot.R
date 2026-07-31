#' @title Volcano plot
#'
#' @description
#' Create a volcano plot from a `SansSouci4DT` object. This volcano plot is made
#' with `ggplot2` and indicates the number of selected proteins with the given
#' thresholds, as well as the minimal number of true positives.
#'
#' @param x A `SansSouci4DT` object. Must be calibrated using the
#' function `sanssouci::fit`.
#' @param pval_thr A `numeric(1)` which is the p-value threshold under which
#' proteins are selected.
#' @param qval_thr A `numeric(1)` which is the q-value threshold under which
#' proteins are selected.
#' @param logfc_thr A `numeric(1)` which is the absolute log(fold change)
#' threshold above which proteins are selected.
#' @param conditions A `character(2)` containing the names of both conditions.
#' @param pal A `character(2)` containing the colors used for non-differential
#' and differential proteins.
#' @param interactive A `logical(1)` defining if the volcano plot should be
#' interactive and use plotly (TRUE) or not and use ggplot2 (FALSE).
#' @param ... unused
#'
#' @returns If `interactive` is TRUE, a `plotly` object providing a volcano
#' plot. Else, a `ggplot2` object providing a volcano plot.
#'
#' @import ggplot2
#' @import sanssouci
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom stats p.adjust
#'
#' @author Manon Gaudin, Nicolas Enjalbert Courrech
#'
#' @exportS3Method sanssouci::volcanoPlot
#'
#' @examples
#' library(sanssouci)
#' n <- 20
#' p <- 32
#' X <- as.data.frame(matrix(rnorm(n * p), ncol = n))
#' colnames(X) <- paste0("S", seq_len(n))
#' X[, "names"] <- paste0("Prot", seq_len(p))
#' group <- rep(c(1, 0), length.out = n)
#' obj <- QFeatures::readQFeatures(X,
#'   quantCols = seq_len(n),
#'   fnames = "names",
#'   name = "datatest"
#' )
#'
#' coldata <- data.frame(
#'   quantCols = paste0("S", seq_len(n)),
#'   Condition = group,
#'   Bio.Rep = as.character(seq_len(n))
#' )
#' rownames(coldata) <- coldata$quantCols
#' SummarizedExperiment::colData(obj) <- coldata
#'
#' ss_obj <- SansSouci4DT(obj)
#' ss_obj <- sanssouci::fit(ss_obj, alpha = 0.5, B = 100)
#'
#' volcanoPlot(ss_obj, pval_thr = 0.5, logfc_thr = 0.05)
#' \dontrun{
#' volcanoPlot(ss_obj, pval_thr = 0.5, logfc_thr = 0.05, interactive = TRUE)
#' }
#'
volcanoPlot.SansSouci4DT <- function(x,
                                     pval_thr = 1,
                                     qval_thr = 1,
                                     logfc_thr = 0,
                                     conditions = NULL,
                                     pal = NULL,
                                     interactive = FALSE, ...) {
  if (missing(x)) {
    stop("'x' is required.")
  }
  sanssouci_object <- x
  if (is.null(sanssouci::pValues(sanssouci_object))) {
    stop("'sanssouci_object' has no p-values associated.
         ('pValues(sanssouci_object)' returns NULL)")
  }
  if (is.null(sanssouci::foldChanges(sanssouci_object))) {
    stop("'sanssouci_object' has no fold-changes associated.
         ('foldChanges(sanssouci_object)' returns NULL)")
  }
  if (length(sanssouci::pValues(sanssouci_object)) !=
        length(sanssouci::foldChanges(sanssouci_object))) {
    stop("p-values and fold-changes associated to 'sanssouci_object'
         must be of equal length.
    ('pValues(sanssouci_object)' and 'foldChanges(sanssouci_object)'
         are not of equal length)")
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
  if (!is.numeric(qval_thr)) {
    stop("'qval_thr' must be numeric.")
  }
  if (length(qval_thr) != 1) {
    stop("'qval_thr' must be of length 1.")
  }
  if ((qval_thr < 0) || (qval_thr > 1)) {
    stop("'qval_thr' must be between 0 and 1.")
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
  if (!is.logical(interactive)) {
    stop("'interactive' must be logical.")
  }
  if (length(interactive) != 1) {
    stop("'interactive' must be of length 1.")
  }

  if ((pval_thr < 1) && (qval_thr < 1)) {
    warning("Filtering both on p-values and BH-adjusted p-values")
  }

  if (is.null(pal)) {
    pal <- c(Out = "gray", In = "orange", "white")
  } else if (length(pal) != 2) {
    warning("'pal' must be of length 2. Colors set to default.")
    pal <- c(Out = "gray", In = "orange", "white")
  } else {
    pal <- c(Out = pal[1], In = pal[2], "white")
  }

  foldchanges <- as.vector(sanssouci::foldChanges(sanssouci_object))
  p.values <- as.vector(sanssouci::pValues(sanssouci_object))
  names(p.values) <- rownames(sanssouci_object$input$Y)

  p <- volcanoPlot(
    x = foldchanges,
    p_value = p.values,
    thr = sanssouci::thresholds(sanssouci_object),
    p_value_bound = p.values,
    p = pval_thr,
    q = qval_thr,
    r = logfc_thr,
    pch = 20,
    cex = c(2, 2),
    col = pal,
    feature_label = "protein",
    add_signed_selections = FALSE
  )
  if (!is.null(conditions) && length(conditions) == 2) {
    p@labels$title <- paste0(conditions[1], "_vs_", conditions[2],
      " - ", p@labels$title,
      sep = ""
    )
  } else if (!is.null(conditions) && length(conditions) != 2) {
    warning("'conditions' must be of length 2. No title added.")
  }

  if (interactive) {
    p <- plotly::ggplotly(p + ggplot2::labs(y = "p-value (-log10 scale)"),
      tooltip = "text"
    ) |>
      plotly::layout(
        margin = list(t = 60, b = 70)
      )
  }
  return(p)
}
