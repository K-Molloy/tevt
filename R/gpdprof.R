#' GPD Profile Log Likelihood for \eqn{\xi}
#'
#' Creates \eqn{\xi} profile log liklihood
#' @param z dataset
#' @param xlow lower \eqn{\xi} value
#' @param xup  upper \eqn{\xi} value
#' @param conf Confidence size \code{default} is set to \code{0.95}
#' @param nint Number of intervals \code{default} is set to \code{100}
#' @examples
#' data(danish)
#' out = gpd.fit(danish, 10)
#' gpd.profxi(out, -0.1, 0.2)
#' @export
gpd.profxi = function(z, xlow, xup, conf = 0.95, nint = 100){
    cat("If routine fails, try changing plotting interval", fill = TRUE)
    xdat <- z$exceedances ; u <- z$threshold
    v <- numeric(nint)
    x <- seq(xup, xlow, length = nint)
    sol <- z$par.ests[1]
    gpd.plikxi <- function(a) {
      # calculates profile log lik
      if(abs(xi) < 10^(-4)) l <- length(xdat) * log(a) + sum(xdat - u)/a
      else {
        y <- (xdat - u)/a
        y <- 1 + xi * y
        if(any(y <= 0) || a <= 0)
          l <- 10^6
        else l <- length(xdat) * log(a) + sum(log(y)) * (1/xi + 1)
      }
      l
    }
    for(i in 1:nint) {
      xi <- x[i]
      opt <- optim(sol, gpd.plikxi, method = "BFGS")
      sol <- opt$par ; v[i] <- opt$value
    }
    plot(x,  - v, type = "l", xlab = "Shape Parameter", ylab =
           "Profile Log-likelihood")
    ma <-  - z$nllh
    abline(h = ma, lty = 1)
    abline(h = ma - 0.5 * qchisq(conf, 1), lty = 1)
    invisible()
  }



#' GPD Profile Log Likelihood for m year return level
#'
#' Creates m year return level plot for profile log liklihood
#' @param z dataset
#' @param m number of years
#' @param xlow lower xi value
#' @param xup  upper xi value
#' @param npy number of days in year \code{default} is set to \code{365}
#' @param conf Confidence size \code{default} is set to \code{0.95}
#' @param nint Number of intervals \code{default} is set to \code{100}
#' @examples
#' data(danish)
#' out = gpd.fit(danish, 10)
#' gpd.prof(out, m = 10, 55, 75)
#' @export
gpd.prof = function(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100){
    cat("If routine fails, try changing plotting interval", fill = TRUE)
    xdat <- z$exceedances ; u <- z$threshold ; la <- z$rate
    v <- numeric(nint)
    x <- seq(xlow, xup, length = nint)
    m <- m * npy
    sol <- z$par.ests[2]
    gpd.plik <- function(a) {
      # calculates profile neg log lik
      if(m != Inf) sc <- (a * (xp - u))/((m * la)^a - 1) else sc <- (u - xp)/
          a
      if(abs(a) < 10^(-4))
        l <- length(xdat) * log(sc) + sum(xdat - u)/sc
      else {
        y <- (xdat - u)/sc
        y <- 1 + a * y
        if(any(y <= 0) || sc <= 0)
          l <- 10^6
        else l <- length(xdat) * log(sc) + sum(log(y)) * (1/a + 1)
      }
      l
    }
    for(i in 1:nint) {
      xp <- x[i]
      opt <- optim(sol, gpd.plik, method = "BFGS")
      sol <- opt$par ; v[i] <- opt$value
    }
    plot(x,  - v, type = "l", xlab = "Return Level", ylab =
           "Profile Log-likelihood")
    ma <-  - z$nllh
    abline(h = ma)
    abline(h = ma - 0.5 * qchisq(conf, 1))
    invisible()
  }
