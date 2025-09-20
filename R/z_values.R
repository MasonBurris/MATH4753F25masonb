#' Standardize a numeric vector
#'
#'This function takes a numeric vector and returns the standardized
#'values
#' @param x A numeric vector.
#'
#' @returns A numeric vector of z values
#' @importFrom stats sd
#' @export
#'
#' @examples
#' mpg <- c(21, 22, 19, 30, 15)
#' z_values(mpg)
z_values <- function(x){
  (x - mean(x, na.rm = TRUE))/(sd(x, na.rm = TRUE))
}
