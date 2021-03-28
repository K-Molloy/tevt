#' Minimizing the distance between the empirical tail and a theoretical Pareto tail with respect to k.
#'
#' An Implementation of the procedure proposed in Danielsson et al. (2016) for selecting the optimal threshold in extreme value analysis.
#' @param data vector of sample data
#' @param ts size of the upper tail the procedure is applied to. Default is 15 percent of the data
#' @param method should be one of \code{ks} for the "Kolmogorov-Smirnov" distance metric or \code{mad} for the mean absolute deviation (default)
#' @details The procedure proposed in Danielsson et al. (2016) minimizes the distance between the largest upper order statistics of the dataset, i.e. the empirical tail, and the theoretical tail of a Pareto distribution. The parameter of this distribution are estimated using Hill's estimator. Therefor one needs the optimal number of upper order statistics \code{k}. The distance is then minimized with respect to this \code{k}. The optimal number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by saying \code{u} is the \code{n-k0}th upper order statistic. For the distance metric in use one could choose the mean absolute deviation called \code{mad} here, or the maximum absolute deviation, also known as the "Kolmogorov-Smirnov" distance metric (\code{ks}). For more information see references.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index by plugging in \code{k0} into the hill estimator}
#' @references Danielsson, J. and Ergun, L.M. and de Haan, L. and de Vries, C.G. (2016). Tail Index Estimation: Quantile Driven Threshold Selection.
#' @examples
#' data(danish)
#' gpd.md(danish,method="mad")
#' gpd.md(danish,method="ks")
#' @export
gpd.md = function(data,ts=0.15,method="mad"){
    xstat=sort(data,decreasing=TRUE)
    n=length(data)
    T=floor(n*ts)
    i=1:(n-1)
    h=(cumsum(log(xstat[i]))/i)-log(xstat[i+1])
    xstat=sort(data)
    A=matrix(ncol=T-1,nrow=T-1)
    for (k in 1:(T-1)){
      for (j in 1:(T-1)){
        A[k,j]=abs( (((k/j)*xstat[n-k+1]^(1/h[k]))^h[k]) - xstat[n-j] )
      }
    }
    if (method=="mad"){
      M=rowMeans(A)
      kstar=which.min(M)
      u=rev(xstat)[kstar]
      list=list(k0=kstar,threshold=u,tail.index=1/h[kstar])
      return(list)
    }
    if (method=="ks"){
      rowMax <- function (rowData) {
        apply(rowData, MARGIN=c(1), max)
      }
      M=rowMax(A)
      kstar=which.min(M)
      u=rev(xstat)[kstar]
      list=list(k0=kstar,threshold=u,tail.index=1/h[kstar])
      return(list)
    }
    if (method!="mad" && method!="ks"){
      warning("method should be one of 'mad' or 'ks'")
    }
  }
