#' recommended value for B
#'
#' @description
#' Provides a recommended value for B, the number of permutation. It gives the
#' minimum value between 1000 and the maximal number of possible permutation for
#' the design in group. The maximal number of possible permutation is
#' \deqn{\binom{n_1 + n_0}{n_1}}, with \(n_1\) and \(n_0\) respectively the
#' number of 1 and 0 in the vector.
#'
#'
#' @param group a vector of '1' and '0
#'
#' @returns a numeric, the recommended value of B
#' @export
#'
#' @examples
#' choose_B(rep(c("1", "0"), times = c(10, 4)))
choose_B <- function(group){
  categCheck(group, length(group))
  n <- length(group)
  n1 <- sum(as.character(group) == "1")
  return(min(1000, choose(n, n1)))
}
