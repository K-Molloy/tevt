#' Generalised Pareto Distribution (GPD)

#' Density, distribution function, quantile function and random number generation
#' for the Generalized Pareto distribution with location, scale, and shape
#' parameters.
#'
#' @name gpd
#' @rdname gpd
#' @param x Vector of observations.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale, and shape parameters. Can be vectors, but
#' the lengths must be appropriate.
#' @param log.d Logical; if TRUE, the log density is returned.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @references Brian Bader, Jun Yan. "eva: Extreme Value Analysis with Goodness-of-Fit Testing." R package version (2016)
#' @examples
#' dgpd(2:4, 1, 0.5, 0.01)
#' dgpd(2, -2:1, 0.5, 0.01)
#' pgpd(2:4, 1, 0.5, 0.01)
#' qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.01)
#' rgpd(6, 1, 0.5, 0.01)
#'
#' ## Generate sample with linear trend in location parameter
#' rgpd(6, 1:6, 0.5, 0.01)
#'
#' ## Generate sample with linear trend in location and scale parameter
#' rgpd(6, 1:6, seq(0.5, 3, 0.5), 0.01)
#'
#' p = (1:9)/10
#' pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#'
#' ## Incorrect syntax (parameter vectors are of different lengths other than 1)
#' # rgpd(1, 1:8, 1:5, 0)
#'
#' ## Also incorrect syntax
#' # rgpd(10, 1:8, 1, 0.01)
#'
#' @details The Generalized Pareto distribution function is given (Pickands, 1975)
#' by \deqn{H(y) = 1 - \Big[1 + \frac{\xi (y - \mu)}{\sigma}\Big]^{-1/\xi}} defined
#' on \eqn{\{y : y > 0, (1 + \xi (y - \mu) / \sigma) > 0 \}}, with location \eqn{\mu},
#' scale \eqn{\sigma > 0}, and shape parameter \eqn{\xi}.
#'
#' @references Pickands III, J. (1975). Statistical inference using extreme order statistics. Annals of Statistics, 119-131.
NULL

#' @rdname gpd
#' @export
dgpd = function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE) {
  # density of the gpd

  # if any scale negative, throw error
  if(min(scale) <= 0)
    stop("Invalid scale")

  # check shape
  cond1 = (length(x) > 1) &
    (((length(loc) != length(x)) & (length(loc) != 1)) |
       ((length(scale) != length(x)) & (length(scale) != 1)) |
       ((length(shape) != length(x)) & (length(shape) != 1)))
  # check shape again
  cond2 = (length(x) == 1) &
    (length(unique(c(length(x), length(loc), length(scale), length(shape)))) > 2)

  # if shape incorrect, throw error
  if(cond1 | cond2)
    stop("Invalid parameter length!")

  # if 1-d, reshape
  if(length(shape) == 1)
    shape = rep(shape, max(length(x), length(loc), length(scale)))

  # calculate density
  below.support = x < loc
  x = pmax(x, loc)
  x = ifelse(shape >= 0, x, pmin(x, (loc - scale/shape)))
  w = (x - loc) / scale
  log.density = -log(scale) - ifelse(shape == 0, w, ((1/shape) + 1) * log1p(w * shape))
  log.density[is.nan(log.density) | is.infinite(log.density) | below.support] = -Inf

  # if exponential is wanted
  if(!log.d)
    log.density = exp(log.density)

  # return density
  log.density
}

#' @rdname gpd
#' @export
pgpd = function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  # distribution function of gpd

  # break if any scale are negative
  if(min(scale) <= 0)
    stop("Invalid scale")

  # do this bit
  cond1 = (length(q) > 1) &
    (((length(loc) != length(q)) & (length(loc) != 1)) |
       ((length(scale) != length(q)) & (length(scale) != 1)) |
       ((length(shape) != length(q)) & (length(shape) != 1)))

  # do this bit too
  cond2 = (length(q) == 1) &
    (length(unique(c(length(q), length(loc), length(scale), length(shape)))) > 2)

  # if cond1 or cond2 are false, some input is wrong length
  if(cond1 | cond2)
    stop("Invalid parameter length!")

  # if the shape is 1-d, reshape it as 2d
  if(length(shape) == 1)
    shape = rep(shape, max(length(q), length(loc), length(scale)))

  # calculate distribution
  q = pmax(q, loc)
  q = ifelse(shape >= 0, q, pmin(q, (loc - scale/shape)))
  w = (q - loc) / scale
  p = ifelse(shape == 0, 1 - exp(-w), 1 - exp((-1/shape)*log1p(w*shape)))

  # if calculating tail
  if(!lower.tail)
    p = 1 - p
  if(log.p)
    p = log(p)

  # return distribution function
  p
}



#' @rdname gpd
#' @export
qgpd = function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  # quantile function for gpd

  # if probabilities are to be given as log(p)
  if(log.p)
    p = exp(p)
  # if probabilities are not in bound [0,1] throw error
  if((min(p, na.rm = TRUE) < 0) || (max(p, na.rm = TRUE) > 1))
    stop("`p' must contain probabilities in (0,1)")

  # if any scale negative, throw error
  if(min(scale) <= 0)
    stop("Invalid scale")

  # if shape is wrong
  cond1 = (length(p) > 1) &
    (((length(loc) != length(p)) & (length(loc) != 1)) |
       ((length(scale) != length(p)) & (length(scale) != 1)) |
       ((length(shape) != length(p)) & (length(shape) != 1)))
  cond2 = (length(p) == 1) &
    (length(unique(c(length(p), length(loc), length(scale), length(shape)))) > 2)
  # if shape is wrong throw error
  if(cond1 | cond2)
    stop("Invalid parameter length!")

  # if probabilities are P[X\leqx], otherwise P[X>x]
  if(lower.tail)
    p = 1 - p
  # if 1-d, reshape
  if(length(shape) == 1)
    shape = rep(shape, max(length(p), length(loc), length(scale)))

  # if gpd shape = 0 (vectorisation is required)
  ifelse(shape == 0, loc - scale * log(p), loc + scale * expm1(-shape * log(p)) / shape)
}


#' @rdname gpd
#' @export
rgpd = function(n, loc = 0, scale = 1, shape = 0) {
  # random generation for gpd

  # if any scale negative, throw error
  if(min(scale) <= 0)
    stop("Invalid scale")

  # check shape
  cond1 = (n > 1) &
    (((length(loc) != n) & (length(loc) != 1)) |
       ((length(scale) != n) & (length(scale) != 1)) |
       ((length(shape) != n) & (length(shape) != 1)))
  cond2 = (n == 1) &
    (length(unique(c(n, length(loc), length(scale), length(shape)))) > 2)
  # if shape incorrect, throw error
  if(cond1 | cond2)
    stop("Invalid parameter length!")

  # run quantilefunction and return
  qgpd(stats::runif(n), loc, scale, shape)
}

#' @rdname gpd
#' @export
experimental.dgpd = function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE) {
  # density of the gpd

  # if any scale negative, throw error
  if(min(scale) <= 0)
    stop("Invalid scale")

  # check shape
  cond1 = (length(x) > 1) &
    (((length(loc) != length(x)) & (length(loc) != 1)) |
       ((length(scale) != length(x)) & (length(scale) != 1)) |
       ((length(shape) != length(x)) & (length(shape) != 1)))
  # check shape again
  cond2 = (length(x) == 1) &
    (length(unique(c(length(x), length(loc), length(scale), length(shape)))) > 2)

  # if shape incorrect, throw error
  if(cond1 | cond2)
    stop("Invalid parameter length!")

  # if 1-d, reshape
  if(length(shape) == 1)
    shape = rep(shape, max(length(x), length(loc), length(scale)))


  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  n = length(x)

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  x = rep(x, length.out = n)
  u = rep(loc, length.out = n)
  sigmau = rep(scale, length.out = n)
  xi = rep(shape, length.out = n)

  d = x # will pass through NA/NaN as entered

  yu = (x - u)/sigmau # used when shape is zero
  syu = 1 + xi*yu     # used when shape non-zero

  # check for x values in range
  yind = ((yu > 0) & (syu > 0))

  d[which(!yind)] = log(0) # zero density is default

  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)

  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    d[whichexp] = -log(sigmau[whichexp]) - yu[whichexp]
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    d[whichxi] = -log(sigmau[whichxi]) - (1/xi[whichxi] + 1) * log(syu[whichxi])
  }

  # return log if required
  if (!log.d) d = exp(d)

  d
}



