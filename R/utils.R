## Find threshold for given number of exceedances
findthresh <- function(data, ne) {
  data <- rev(sort(as.numeric(data)))
  data[ne] - min(min(abs(diff(data))[abs(diff(data)) > 0]), 1e-6)
}


## Function to help deal with design matrix
adjScale = function (x)
{
  truemeans <- as.numeric(colMeans(x))
  truevars <- as.numeric(apply(x, 2, sd))
  adjmeans <- ifelse(truevars == 0, 0, truemeans)
  adjvars <- ifelse(truevars == 0, truemeans, truevars)
  if (ncol(x) == 1)
    adjmeans <- 0
  x <- t((t(x) - adjmeans)/adjvars)
  out <- list(x, truemeans, truevars, adjmeans, adjvars)
  names(out) <- c("mat", "truemeans", "truevars", "adjmeans",
                  "adjvars")
  out
}


## Returns expected inverse fisher information matrix
gpdFisher = function (n, theta)
{
  scale <- theta[1]
  shape <- theta[2]
  one <- (2 * (1 + shape) * scale^2)/n
  two <- (1 + shape)^2/n
  cov <- -((1 + shape) * scale)/n
  varcov <- matrix(c(one, cov, cov, two), 2)
  varcov
}


## Stopping criterion for sequential testing
pSeqStop = function (p)
{
  m <- length(p)
  int <- seq(1, m, 1)
  pFDR <- cumsum(-log(1 - p[int]))/int
  pFWER <- rev(exp(cumsum(rev(log(p[int]))/rev(int)))) * (m/int)
  out <- cbind.data.frame(pFWER, pFDR, p)
  colnames(out) <- c("StrongStop", "ForwardStop", "UnAdjusted")
  out
}
