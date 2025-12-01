#' ntickets
#'
#' Package for overbooking problem.
#'
#' @param N integer, number of seats.
#' @param gamma probability tolerance for true overbooking.
#' @param p actual attendees probability.
#' @param make_plots plot discrete and continuous objectives.
#'
#' @returns named list containing nd, nc, N, p and gamma - where nd is
#' calculated using the discrete distribution and nc is the same calculated with
#' normal approximation.
#' @export
#'
#' @examples
#' ntickets(400, 0.02, 0.95)
ntickets <- function(N, gamma, p, make_plots = TRUE) {
  fd <- function(n) pbinom(N, n, p) - (1 - gamma)
  fc <- function(n) {
    z <- (N + 0.5 - n * p) / sqrt(n * p * (1 - p))
    pnorm(z) - (1 -gamma)
  }

  grid <- seq(max(N, floor(N/p - 15)), ceiling(N/p + 15), by = 1)

  vals_d <- fd(grid)
  nd <- if (any(vals_d >= 0)) max(grid[vals_d >= 0]) else NA_integer_

  hi <- max(grid)
  nc_root <- try(uniroot(function(x) fc(x), interval = c(N, hi))$root, silent = TRUE)

  if(inherits(nc_root, "try-error")) {
    vals_c <- fc(grid)
    nc <- grid[which.min(abs(vals_c))]
  } else {
    nc <- floor(nc_root)
  }

  if (make_plots) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(1,2), mar = c(4,4,2,1))

    # Discrete plot
    plot(grid, abs(vals_d), type = "b", pch = 21, bg = "blue",
         xlab = "n (tickets sold)", ylab = "Objective",
         main = "Discrete Objective")
    abline(h = 0, lty = 2)
    if(!is.na(nd)) points(nd, fd(nd), pch = 21, bg = "red", cex = 1.2)

    # Continuous Plot
    vals_c <- fc(grid)
    plot(grid, abs(vals_c), type = "l", lwd = 2,
         xlab = "n (tickets sold)", ylab = "Objective",
         main = "Normal approx objective")
    abline(h = 0, lty = 2)
    if(!is.na(nc)) points(nd, fd(nd), pch = 21, bg = "blue", cex = 1.2)
  }

  list(nd = nd, nc_root = nc_root, N = N, p = p, gamma = gamma)
}







