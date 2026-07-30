#' @title Limma test for rows of a matrix
#'
#' @param X A \code{m x n} `numeric matrix` whose rows correspond to variables
#'   and columns to observations.
#' @param categ Either a `numeric vector` of \code{n} categories in \eqn{0, 1}
#'   for the observations, or a \code{n x B} matrix stacking \code{B} such vectors
#'   (typically permutations of an original vector of size \code{n}).
#' @param alternative DEPRECATED  A `character string` specifying the
#'   alternative hypothesis "two.sided".
#'
#' @returns A `list` containing the following components:
#' \describe{
#'   \item{p.value}{the p-values for the tests}
#'   \item{estimate}{the mean difference between groups}}
#'   Each of these elements is a matrix of size \code{m x B}, coerced to a
#'   vector of length \code{m} if \code{B=1}
#'
#' @author Nicolas Enjalbert Courrech
#'
#' @export
#'
#' @examples
#' m <- 300
#' n <- 38
#' mat <- matrix(rnorm(m * n), ncol = n)
#' categ <- rep(c(0, 1), times = c(27, n - 27))
#' res <- rowLimmaTest(mat, categ, alternative = "two.sided")
#' 
#' #If categ is a matrix
#' categ_mat <- replicate(10, sample(categ))
#' res <- rowLimmaTest(mat, categ_mat, alternative = "two.sided")
rowLimmaTest <- function(X, categ, alternative = c("two.sided")) {
  categ <- as.matrix(categ)
  B <- ncol(categ)
  p.value <- matrix(NA, nrow = nrow(X), ncol = B)
  logFC <- matrix(NA, nrow = nrow(X), ncol = B)
  for (b in seq(B)) {
    categ_b <- categ[, b]
    res <- limmaAnalyses(Y = X, groups = categ_b)
    p.value[, b] <- res$p.values
    logFC[, b] <- res$logFC
  }
  if (B == 1) {
    p.value <- as.vector(p.value)
    logFC <- as.vector(logFC)
  }
  return(list(p.value = p.value, estimate = logFC))
}

#' @title Wrapper of limma Differential Expression Analysis pipeline
#'
#' @param Y A matrix (p*n) with variables to test in rows and samples in column.
#' @param groups A binary vector (of size n) of two conditions samples.
#'
#' @returns A list with 'p.values' a vector of size p of raw pvalues
#'  and 'logFC' a vector of size p of log Fold Change, both from Limma pipeline.
#'
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom stats model.matrix
#'
#' @author Nicolas Enjalbert Courrech
#'
#' @keywords internal
#' @export
#'
#' @examples
#' Y <- matrix(rnorm(50), ncol = 5)
#' groups <- c(1, 1, 0, 0, 0)
#' res <- limmaAnalyses(Y, groups)
#'
limmaAnalyses <- function(Y, groups) {
  categCheck(groups, ncol(Y))
  groups <- as.factor(groups)
  # Create design matrix with no intercepts
  design.matrix <- stats::model.matrix(~ 0 + groups)
  # Fit a linear model
  res_lm <- limma::lmFit(Y, design.matrix)
  # Define contrast : here onlye group1 vs group0
  contr <- limma::makeContrasts("groups1 - groups0",
    levels = colnames(design.matrix)
  )
  # fit with contrast
  res_fit <- limma::contrasts.fit(res_lm, contr)
  # make test
  res_eb <- limma::eBayes(res_fit)

  return(list(p.values = res_eb$p.value, logFC = res_eb$coefficients))
}
