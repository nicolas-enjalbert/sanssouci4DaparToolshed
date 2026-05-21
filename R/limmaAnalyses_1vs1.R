#' Wrapper of limma DEA pipeline
#'
#' @param Y a matrix (p*n) with variables to test in rows and samples in column
#' @param groups a binary vector (of size n) of two conditions samples
#'
#' @returns a list with 'p.values' a vector of size p of raw pvalues
#'  and 'logFC' a vector of size p of log Fold Change, both from Limma pipeline
#'
#'  @importFrom limma lmFit
#'  @importFrom limma makeContrasts
#'  @importFrom limma contrasts.fit
#'  @importFrom limma eBayes
#'  @importFrom stats model.matrix
#'
#'
#' @keywords internal
#' @export
#'
#' @examples
#' Y <- matrix(rnorm(50), ncol = 5)
#' groups = c(1,1,0,0,0)
#' res <- limmaAnalyses_1vs1(Y, groups)
limmaAnalyses_1vs1 <- function(Y, groups){
  groups <- as.factor(groups)
  # Create design matrix with no intercepts
  design.matrix <- stats::model.matrix(~ 0 + groups)
  # Fit a linear model
  res_lm <- limma::lmFit(Y, design.matrix)
  # Define contrast : here onlye group1 vs group0
  contr <- limma::makeContrasts(groups1 - groups0,
                                levels = colnames(design.matrix))
  # fit with contrast
  res_fit <- limma::contrasts.fit(res_lm, contr)
  # make test
  res_eb <- limma::eBayes(res_fit)
  #extract
  # TT <- limma::topTable(res_eb, sort.by = "none", number = Inf)

  return(list(p.values = res_eb$p.value, logFC = res_eb$coefficients))
}
