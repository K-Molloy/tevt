#' Sum Plot
#'
#' An implementation of the so called sum plot proposed in de Sousa & Michailidis (2004)
#' @param data vector of sample data
#' @param kmin gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}.
#' @details The sum plot is based on the plot \code{(k,S_k)} with \code{S_k:=k*gamma_k} where \code{gamma_k} denotes the Hill estimator. So the sum plot and the Hill plot are statistically equivalent. The sum plot should be approximately linear for the \code{k}-values where \code{gamma_k=gamma}. So the linear part of the graph can be used as an estimator of the (inverse) tail index. The sum plot leads to the estimation of the slope while the classical Hill plot leads to estimation of the intercept. The optimal number of order statistics, also known as the threshold, can then be derived as the value \code{k} where the plot differs from a straight line with slope \code{gamma}. See references for more information.
#' @return
#' The plot shows the values of \code{S_k=k*gamma_k} for different \code{k}. See references for more information.
#' @references De Sousa, Bruno and Michailidis, George (2004). A diagnostic plot for estimating the tail index of a distribution. \emph{Journal of Computational and Graphical Statistics} \bold{13(4)}, 1--22.
#' @examples
#' data(danish)
#' gpd.plot.sum(danish)
#' @export
gpd.plot.sum=function (data, kmin = 5){

  # get size of data
  n = length(data)

  i = 1:(n - 1)
  # sort it
  xstat = sort(data, decreasing = TRUE)
  # calculate h
  h = (cumsum(log(xstat[i]))/i) - log(xstat[i + 1])
  # plot it
  plot(i[kmin:(n - 1)], i[kmin:(n-1)]/h[kmin:(n - 1)], type = "l", xlab = "Order Statistics",
       ylab = "Tail Index")
}
