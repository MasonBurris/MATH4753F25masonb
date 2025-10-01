#' Birthday Problem
#'
#' Computes the probability of two people in a group of size n
#' sharing the same birthday.
#'
#' @param n Integer of class size.
#'
#' @returns A numeric vector of probabilities.
#' @export
#'
#' @examples
#' birthday(23)
#' birthday(1:20)
birthday <- function(n) {

  1 - exp(lfactorial(n) + lchoose(365, n) - n * log(365))

}
