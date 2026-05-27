#' Check design vector
#'
#' @param categ a vector of design. Should consist only of '0' and '1'. If not return an error.
#' @param n expected length of categ
#'
#' @keywords internal
#'
#' @returns NULL
categCheck <- function(categ, n) {
  name <- as.character(substitute(categ))
  if (length(categ) != n) {
    stop(name, " should be of length ", n, ", not ", length(categ))
  }
  categ <- as.factor(categ)
  cats <- levels(categ)
  if (!identical(cats, c("0", "1")) & length(cats) <= 2) {
    stop("'", name, "' should consist only of '0' and '1' or distinct continuous values.")
  }
}
