#' @title Volcano plot
#'
#' @description
#' Create a volcano plot from a `Sanssouci` object. This volcano plot is made
#' with `ggplot2` and indicates the number of selected proteins with the given
#' thresholds, as well as the minimal number of true positives.
#'
#' @param sanssouci_object A `SansSouci` object. Must be calibrated using the
#' function `fit`.
#' @param pval_thr A `numeric(1)` which is the p-value threshold under which
#' proteins are selected.
#' @param qval_thr A `numeric(1)` which is the q-value threshold under which
#' proteins are selected.
#' @param logfc_thr A `numeric(1)` which is the absolute log(fold change)
#' threshold above which proteins are selected.
#' @param conditions A `character(2)` containing the names of both conditions.
#' @param pal A `character(2)` containing the colors used for differential and
#' non differential proteins.
#' @param interactive A `logical(1)` defining if the volcano plot should be
#' interactive and use plotly (TRUE) or not and use ggplot2 (FALSE).
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
#' VolcanoPlot_ss4DT(
#'   sanssouci_object = ss_obj_fit,
#'   pval_thr = 0.5,
#'   logfc_thr = 0.05
#' )
#'
VolcanoPlot_ss4DT <- function(sanssouci_object,
                              pval_thr = 1,
                              qval_thr = 1,
                              logfc_thr = 0,
                              conditions = NULL,
                              pal = NULL,
                              interactive = FALSE) {
  if (missing(sanssouci_object)) {
    stop("'sanssouci_object' is required.")
  }
  if (!inherits(sanssouci_object, "SansSouci")) {
    stop("'sanssouci_object' must be an object of class SansSouci.")
  }
  if (is.null(sanssouci::pValues(sanssouci_object))) {
    stop("'sanssouci_object' has no p-values associated.
         ('pValues(sanssouci_object)' returns NULL)")
  }
  if (is.null(sanssouci::foldChanges(sanssouci_object))) {
    stop("'sanssouci_object' has no fold-changes associated.
         ('foldChanges(sanssouci_object)' returns NULL)")
  }
  if (length(sanssouci::pValues(sanssouci_object)) != length(sanssouci::foldChanges(sanssouci_object))) {
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

  if (is.null(conditions)) {
    title <- ""
  } else if (!is.null(conditions) && length(conditions) != 2) {
    warning("'conditions' must be of length 2. No title added.")
    title <- ""
  } else {
    title <- paste0(conditions[1], "_vs_", conditions[2])
  }

  if (is.null(pal)) {
    pal <- list(In = "orange", Out = "gray")
  } else if (length(pal) != 2) {
    warning("'pal' must be of length 2. Colors set to default.")
    pal <- list(In = "orange", Out = "gray")
  } else {
    pal <- list(In = pal[1], Out = pal[2])
  }

  pval <- as.vector(sanssouci::pValues(sanssouci_object))
  logFC <- as.vector(sanssouci::foldChanges(sanssouci_object))

  logpval <- -log10(pval)

  adjp <- p.adjust(pval, method = "BH") ## adjusted p-values
  pval_sel <- which((adjp <= qval_thr) & ## selected by q-value
                      (pval <= pval_thr)) ##          or p-value
  log_pval_thr <- Inf
  if (length(pval_sel) > 0) {
    log_pval_thr <- min(logpval[pval_sel]) ## threshold on the log(p-value) scale
  }

  ## protein selections
  sel1 <- which(logpval >= log_pval_thr & logFC >= logfc_thr)
  sel2 <- which(logpval >= log_pval_thr & logFC <= -logfc_thr)
  sel12 <- sort(union(sel1, sel2))

  ## post hoc bounds in selections
  pval_post_hoc <- sanssouci::pValues(sanssouci_object)
  thr_post_hoc <- sanssouci::thresholds(sanssouci_object)

  # n1 <- length(sel1)
  # FP1 <- maxFP(pval_post_hoc[sel1], thr = thr_post_hoc)
  # TP1 <- n1 - FP1
  # FDP1 <- round(FP1 / max(n1, 1), 2)
  #
  # n2 <- length(sel2)
  # FP2 <- maxFP(pval_post_hoc[sel2], thr = thr_post_hoc)
  # TP2 <- n2 - FP2
  # FDP2 <- round(FP2 / max(n2, 1), 2)

  n12 <- length(sel12)
  FP12 <- sanssouci::maxFP(pval_post_hoc[sel12], thr = thr_post_hoc)
  TP12 <- n12 - FP12
  FDP12 <- round(FP12 / max(n12, 1), 2)


  names_prot <- rownames(sanssouci_object$input$Y)
  if (is.null(names_prot)) {
    names_prot <- seq_along(logFC)
  }

  df <- data.frame(pval = logpval, logFC = logFC, protein_name = names_prot)

  # Significant / non-significant groups
  df$isDiff <- ifelse(df$pval >= log_pval_thr & abs(df$logFC) >= logfc_thr,
                      "In",
                      "Out"
  )

  if (title != "") {
    title <- paste0(
      title, " - ", n12, " proteins selected \n",
      "At least ", TP12, " true positives (FDP <= ", FDP12, ")"
    )
  } else {
    title <- paste0(
      n12, " proteins selected \n",
      "At least ", TP12, " true positives (FDP <= ", FDP12, ")"
    )
  }

  p <- suppressWarnings(ggplot2::ggplot(
    df,
    ggplot2::aes(x = logFC, y = pval, color = isDiff)
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(pval_thr),
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logfc_thr, logfc_thr),
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "solid",
      color = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(text = paste0(
        "Protein: ", protein_name, "<br>",
        "logFC: ", logFC, "<br>",
        "-log10(pval): ", pval
      )),
      alpha = 0.8,
      size = 2
    ) +
    ggplot2::scale_color_manual(values = c(In = pal$In, Out = pal$Out)) +
    ggplot2::labs(
      title = title,
      x = "logFC",
      y = "-log10(pValue)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      legend.position = "none"
    ))

  if (interactive) {
    p <- plotly::ggplotly(p, tooltip = "text") |>
      plotly::layout(
        title = list(
          text = title,
          font = list(size = 18)
        ),
        margin = list(t = 60, b = 70)
      )
  }

  return(p)
}
