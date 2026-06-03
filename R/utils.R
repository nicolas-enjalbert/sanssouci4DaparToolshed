#' @title Check design vector
#'
#' @description
#' Check whether the design vector is valid. 
#'
#' @param categ A vector of design. Should consist only of '0' and '1'.
#' If not return an error.
#' @param n A `numeric(1)`, expected length of categ.
#'
#' @keywords internal
#'
#' @returns NULL if the design vector is valid, or an error if the design 
#' vector is invalid.
#'
#' @author Nicolas Enjalbert Courrech
#'
#' @export
#'
#' @examples
#' group <- rep(c("1", "0"), times = c(10, 4))
#' categCheck(group, length(group))
#'
categCheck <- function(categ, n) {
  name <- as.character(substitute(categ))
  if (length(categ) != n) {
    stop(name, " should be of length ", n, ", not ", length(categ))
  }
  categ <- as.factor(categ)
  cats <- levels(categ)
  if (!identical(cats, c("0", "1")) && length(cats) <= 2) {
    txt <- paste0("'", name, "' should consist only of '0' and '1' or distinct 
                  continuous values.")
    stop(txt)
  }
}
