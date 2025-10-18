#' Plot Normal Curve with Shaded Area
#'
#'This function graphs a normal distribution with mean mu and standard
#'deviation sigma, shades the area from negative infinity to 'a', and returns
#'the calculated probability.
#'
#' @param mu Mean of the normal distribution.
#' @param sigma Standard deviation of the normal distribution.
#' @param a Upper bound for the shaded area P(X <= a).
#'
#' @returns
#' A named list with:
#' \item{mu}{The mean used in the plot.}
#' \item{sigma}{The standard deviation used in the plot.}
#' \item{area}{The probability P(X <= a).}
#'
#' @export
#'
#' @examples
#' myncurve(mu = 10, sigma = 5, a = 6)
myncurve <- function(mu, sigma, a) {
  # draw the normal pdf
  xmin <- mu - 3*sigma
  xmax <- mu + 3*sigma
  curve(dnorm(x, mean = mu, sd = sigma),
        xlim = c(xmin, xmax),
        main = paste("mu =", mu, ", sigma =", sigma),
        xlab = "x", ylab = "Density")

  # shade from -Inf (clamped at xmin) to a (clamped at xmax)
  x_to <- max(min(a, xmax), xmin)
  xcur <- seq(xmin, x_to, length = 1000)
  ycur <- dnorm(xcur, mean = mu, sd = sigma)
  polygon(c(xmin, xcur, x_to), c(0, ycur, 0), col = "Blue")

  # area = P(X <= a)
  area <- pnorm(a, mean = mu, sd = sigma)
  area <- round(area, 4)

  # annotate (place label a bit to the right of mu)
  text(x = mu + 0.8*sigma, y = dnorm(mu, mu, sigma)*0.8,
       paste("Area = ", area, sep = ""))

  # return named list
  list(mu = mu, sigma = sigma, area = area)
}
