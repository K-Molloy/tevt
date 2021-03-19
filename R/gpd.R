## gpd stuff

# Density, distribution function, quantile function and random number generation
# for the Generalized Pareto distribution with location, scale, and shape
# parameters.

# Conforming to standard R distribution function parameters


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
  if((min(p, na.rm = TRUE) <= 0) || (max(p, na.rm = TRUE) >= 1))
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

  # if gpd shape = 0
  if (shape==0){
    q = loc - scale * log(p)
  } else {
    q = loc + scale * expm1(-shape*log(p)) / shape
  }
  # return quantile
  q
  #ifelse(shape == 0, loc - scale * log(p), loc + scale * expm1(-shape * log(p)) / shape)
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
  qgpd(runif(n), loc, scale, shape)
}



