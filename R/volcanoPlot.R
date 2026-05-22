#' @title Volcano plot
#'
#' @description 
#' Create a volcano plot
#'
#' @param pval A `numeric` vector of p-value values. 
#' This vector must have the same length as `logFC`.
#' @param logFC A `numeric` vector of log(fold change) values. 
#' This vector must have the same length as `pval`.
#' @param th_pval A `numeric(1)` which is the p-value threshold under which 
#' proteins are selected.
#' @param th_logfc A `numeric(1)` which is the absolute log(fold change) 
#' threshold above which proteins are selected.
#' @param conditions A `character(2)` containing the names of both conditions.
#' @param pal A `character(2)` containing the colors used for differential and 
#' non differential proteins. 
#'
#' @returns A volcano plot.
#'
#' @import ggplot2
#'
#' @author Manon Gaudin
#'
#' @export
#'
#' @examples
#' #xxx a faire
#' NULL
#' 

VolcanoPlot_ggplot <- function(pval,
                               logFC,
                               th_pval = 0,
                               th_logfc = 0,
                               conditions = NULL,
                               pal = NULL) {
  if (missing(pval)){
    stop("'pval' is required.")
  }
  if (!is.numeric(pval)){
    stop("'pval' must be an object of class numeric.")
  }
  if (missing(logFC)){
    stop("'logFC' is required.")
  }
  if (!is.numeric(logFC)){
    stop("'logFC' must be an object of class numeric.")
  }
  if (length(pval) != length(logFC)){
    stop("'pval' and 'logFC' must be of equal length.")
  }
  if (!is.numeric(th_pval)){
    stop("'th_pval' must be numeric.")
  }
  if (length(th_pval) != 1){
    stop("'th_pval' must be of length 1.")
  }
  if (th_pval < 0){
    stop("'th_pval' must be positive.")
  }
  if (!is.numeric(th_logfc)){
    stop("'th_logfc' must be numeric.")
  }
  if (length(th_logfc) != 1){
    stop("'th_logfc' must be of length 1.")
  }
  if (th_logfc < 0){
    stop("'th_logfc' must be positive.")
  }
  
  if (is.null(conditions)){
    title <- ""
  } else if (!is.null(conditions) & length(conditions) != 2) {
    warning("'conditions' must be of length 2. No title added.")
    title <- ""
  } else {
    title <- paste0(conditions[1], "_vs_", conditions[2])
  }
  
  if (is.null(pal)) {
    pal <- list(In = "orange", Out = "gray")
  } else if (length(pal) != 2) {
    warning("The palette must be a list of two items: In and Out. Set to default.")
    pal <- list(In = "orange", Out = "gray")
  } else {
    pal <- list(In = pal[1], Out = pal[2])
  }
  
  df <- data.frame(pval = pval, logFC = logFC)
  
  # Significant / non-significant groups
  df$isDiff <- ifelse(df$pval >= th_pval & abs(df$logFC) >= th_logfc,
                      "In", 
                      "Out")
  
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = logFC, y = pval, color = isDiff)) +
    ggplot2::geom_hline(yintercept = th_pval,
                        linetype = "dashed",
                        color = "grey") +
    ggplot2::geom_vline(xintercept = c(-th_logfc, th_logfc),
                        linetype = "dashed",
                        color = "grey") +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = "solid",
                        color = "black") +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "solid",
                        color = "black") +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::scale_color_manual(values = c(In = pal$In, Out = pal$Out)) +
    ggplot2::labs(title = title,
                  x = "logFC",
                  y = "-log10(pValue)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 16),
                   legend.position = "none")
  
  return(p)
}
