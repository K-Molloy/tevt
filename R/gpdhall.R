#' A Single Bootstrap Procedure for Choosing the Optimal Sample Fraction
#'
#' An Implementation of the procedure proposed in Hall (1990) for selecting the optimal sample fraction in tail index estimation
#' @param data vector of sample data
#' @param B number of Bootstrap replications
#' @param epsilon gives the amount of the first resampling size \code{n1} by choosing \code{n1 = n^epsilon}. Default is set to \code{epsilon=0.955}
#' @param kaux tuning parameter for the hill estimator
#' @details The Bootstrap procedure simulates the AMSE criterion of the Hill estimator. The unknown theoretical parameter of the inverse tail index \code{gamma} is replaced by a consistent estimation using a tuning parameter \code{kaux} for the Hill estimator. Minimizing this statistic gives a consistent estimator of the sample fraction \code{k/n} with \code{k} the optimal number of upper order statistics. This number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#'@references Hall, P. (1990). Using the Bootstrap to Estimate Mean Squared Error and Select Smoothing Parameter in Nonparametric Problems. \emph{Journal of Multivariate Analysis}, \bold{32}, 177--203.
#' @examples
#' data(danish)
#' gpd.hall(danish)
#' @export
gpd.hall = function(data,B=1000,epsilon=0.955,kaux=2*sqrt(length(data))){
    n=length(data)
    n1=floor(n^epsilon)

    helphill=function (k) {
      xstat = sort(data, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat
    }

    help=helphill(kaux)

    mse=matrix(nrow=B,ncol=n1-1)
    for (l in 1:B){
      i=1:(n1-1)
      x1=sample(data,n1,replace=TRUE)
      x1=sort(x1,decreasing=TRUE)
      h=(cumsum(log(x1[i]))/i)-log(x1[i+1])
      mse[l,]=h-help
    }
    mse=mse^2
    msestar=colMeans(mse)
    k1star=which.min(msestar)

    k0star=floor(k1star*(n/n1)^(2/3))
    u=sort(data,decreasing=TRUE)[k0star]
    ti=1/helphill(k0star)
    list=list(k0=k0star,threshold=u,tail.index=ti)
    list
  }
