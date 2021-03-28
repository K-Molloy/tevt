#' Diagnostic plots for a fit to the Generalized Pareto distribution
#'
#' Allows calling of a 2x2 diagnostic plot, different plots available for
#' stationary vs non-stationary models.
#' The extension from EVA allows for direct calling to gpd.diag and plots
#' produced. Additionally, it exposes the individual plots for custom sizing
#' (especially via .rmd or .rnw)
#'
#' @param z A class object returned from `gpdFit'.
#' @param conf Confidence level used in the return level plot.
#' @param method The method to compute the return level confidence interval - either delta method (default) or
#' profile likelihood. Choosing profile likelihood may be quite slow.
#' @param safe Binary (default FALSE) For when plot margins are too large
#' @examples
#' \dontrun{
#' x <- rgpd(10000, loc = 0.5, scale = 1, shape = 0.1)
#' z <- gpd.fit(x, nextremes = 500)
#' plot(z)
#' }
#' @return For stationary models, provides return level, density, probability, and quantile plots for the GPD exceedances. The
#' overlaid density is the `true' density for the estimated parameters. For nonstationary models, provides
#' residual probability and quantile plots. In addition, nonstationary models provide plots of the residuals vs.
#' the parameter covariates.
#' @details See the reference for details on how return levels are calculated.
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
#' @references Modified from "eva : Extreme Value Analysis with Goodness-of-Fit Testing."
#' @importFrom utils menu
#' @export
gpd.diag <- function(z, conf = 0.95, method = c("delta", "profile"), safe=FALSE) {
  # make sure par returns to what it used to be
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  method <- match.arg(method)
  par(mfcol = c(2, 2))
  scalevec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  shapevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))

  if (safe) par(mar=c(1,1,1,1))
  if(z$stationary){
    gpd.rlPlot(z, conf, method)
    gpd.hist(z)
  }
  if(!z$stationary){
    gpd.resid(z, scalevec, shapevec)
  }
  gpd.pp(z, scalevec, shapevec)
  gpd.qq2(z, scalevec, shapevec)
}

#' Return level Plot
#' @param z GPD Model
#' @param conf Confidence Level
#' @param method Method of Optimisation type
#' @export
gpd.rlPlot <- function(z, conf = 0.95, method = c("delta", "profile")) {
  if(!z$stationary)
    stop("Model must be stationary")
  method <- match.arg(method)
  emp <- z$n / z$npp
  p <- c(seq(from = 0.001, to = emp, by = emp / 100), seq(from = emp, to = 100 * emp, by = emp))
  levels <- matrix(0, length(p), 3)
  for(i in 1:nrow(levels)) {
    y <- gpd.rl(z, p[i], conf = conf, method = method)
    levels[i, 1] <- y$Estimate
    levels[i, 2:3] <- y$CI
  }
  plot((z$n + 1) / ((1:z$n) * z$npp), rev(sort(z$data)), type = "n", log = "x",
       xlab = "Return Period", ylab = "Return Level", xlim = c(min((z$n + 1) / ((1:z$n) * z$npp)), max(p)),
       ylim = c(min(z$data, levels[, 2]), max(z$data, levels[, 3])))
  title("Return Level Plot")
  lines(p, levels[, 1])
  lines(p, levels[, 2], col = 4)
  lines(p, levels[, 3], col = 4)
  points((z$n + 1) / ((1:z$n) * z$npp), rev(sort(z$data)))
}

#' Histogram
#' @param z GPD Model
#' @export
gpd.hist <- function(z) {
  if(!z$stationary)
    stop("Model must be stationary")
  excess <- z$data[z$data > z$threshold]
  h <- hist(excess, plot = FALSE)
  x <- seq(min(excess), max(excess), (max(excess) - min(excess))/1000)
  y <- dgpd(x, loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2])
  hist(excess, freq = FALSE, xlim = c(min(excess), max(excess)), ylim = c(0, max(max(h$density), max(y))),
       xlab = "x", ylab = "Density", main = "Density Plot")
  points(excess, rep(0, length(excess)))
  lines(x, y, col = 4)
}


#' QQ-Probability Plot
#' @param z GPD Model
#' @param scalevec Vector of scale
#' @param shapevec Vector of Shape
#' @export
gpd.qq2 <- function(z, scalevec, shapevec) {
  # allow direct calling to qq2
  if(missing(scalevec)) scalevec = z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  if(missing(shapevec)) shapevec = z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))

  excess <- (z$data - z$threshold)[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  emp <- qgpd(Series, loc = 0, scale = scalevec, shape = shapevec)
  plot(sort(excess), emp, xlab = "Empirical", ylab = "Model",
       xlim = c(min(excess, emp), max(excess, emp)), ylim = c(min(excess, emp), max(excess, emp)))
  if(z$stationary)
    title("Quantile Plot")
  if(!z$stationary)
    title("Residual Quantile Plot")
  abline(0, 1, col = 4)
}

#' PP-Probability Plot
#' @param z GPD Model
#' @param scalevec Vector of scale
#' @param shapevec Vector of Shape
#' @export
gpd.pp <- function(z, scalevec, shapevec) {
  # allow direct calling to pp
  if(missing(scalevec)) scalevec = z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  if(missing(shapevec)) shapevec = z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))

  excess <- (z$data - z$threshold)[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  p <- pgpd(excess, loc = 0, scale = scalevec, shape = shapevec)
  plot(sort(p), Series, xlab = "Empirical", ylab = "Model", xlim = c(0,1), ylim = c(0,1))
  if(z$stationary)
    title("Probability Plot")
  if(!z$stationary)
    title("Residual Probability Plot")
  abline(0, 1, col = 4)
}

#' residuals vs. the covariates Plot
#' @param z GPD Model
#' @param scalevec Vector of scale
#' @param shapevec Vector of Shape
#' @export
gpd.resid <- function(z, scalevec, shapevec) {
  if(z$stationary)
    stop("Model cannot be stationary")
  # allow direct calling to resid
  if(missing(scalevec)) scalevec = z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  if(missing(shapevec)) shapevec = z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))

  excess <- (z$data - z$threshold)[z$data > z$threshold]
  resid <- pgpd(excess, loc = 0, scale = scalevec, shape = shapevec)
  if(z$parnum[1] > 1) {
    for(i in 2:z$parnum[1]) {
      plot(z$covars[[1]][, i], resid, xlab = paste("Scale", colnames(z$covars[[1]])[i], sep = " "), ylab = "Residuals")
      lines(lowess(z$covars[[1]][, i], resid), col = "red")
    }
  }
  if(z$parnum[2] > 1) {
    for(i in 2:z$parnum[2]) {
      plot(z$covars[[2]][, i], resid, xlab = paste("Shape", colnames(z$covars[[2]])[i], sep = " "), ylab = "Residuals")
      lines(lowess(z$covars[[2]][, i], resid), col = "red")
    }
  }
}
