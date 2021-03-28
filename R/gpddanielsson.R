#' A Double Bootstrap Procedure for Choosing the Optimal Sample Fraction
#'
#' An Implementation of the procedure proposed in Danielsson et al. (2001) for selecting the optimal sample fraction in tail index estimation.
#' @param data vector of sample data
#' @param B number of Bootstrap replications
#' @param epsilon gives the amount of the first resampling size \code{n1} by choosing \code{n1 = n^epsilon}. Default is set to \code{epsilon=0.9}
#' @details The Double Bootstrap procedure simulates the AMSE criterion of the Hill estimator using an auxiliary statistic. Minimizing this statistic gives a consistent estimator of the sample fraction \code{k/n} with \code{k} the optimal number of upper order statistics. This number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{second.order.par}{gives an estimation of the second order parameter \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#'@references Danielsson, J. and Haan, L. and Peng, L. and Vries, C.G. (2001). Using a bootstrap method to choose the sample fraction in tail index estimation. \emph{Journal of Multivariate analysis}, \bold{2}, 226-248.
#' @examples
#' data=rexp(100)
#' gpd.danielsson(data, B=200)
#' @export
gpd.danielsson = function(data,B=500,epsilon=0.9){
    n=length(data)
    n1=floor(n^epsilon)
    n2=(n1^2)/n

    Qn1=function (k) {
      xstat = sort(x1, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
      xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat-(2*xihat2^2)
    }
    Qn2=function (k) {
      xstat = sort(x2, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
      xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat-(2*xihat2^2)
    }

    qn1=matrix(nrow=B,ncol=n1-1)
    qn2=matrix(nrow=B,ncol=n2-1)
    for (l in 1:B){
      x1=sample(data,n1,replace=TRUE)
      x2=sample(data,n2,replace=TRUE)
      qn1[l,]=sapply(1:(n1-1),Qn1)
      qn2[l,]=sapply(1:(n2-1),Qn2)
    }
    qn1=qn1^2
    qn2=qn2^2
    qn1star=colMeans(qn1)
    qn2star=colMeans(qn2)
    k1star=which.min(qn1star)
    k2star=which.min(qn2star)

    Exp=(log(n1)-log(k1star))/(log(n1))
    Z=(log(k1star))^2
    N=(2*log(n1)-log(k1star))^2
    k0star=floor(k1star^2/k2star*((Z/N)^Exp))+1
    u=sort(data,decreasing=TRUE)[k0star]
    rho=log(k1star)/(-2*log(n1)+2*log(k1star))

    helphill=function (k) {
      xstat = sort(data, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat
    }

    ti=1/helphill(k0star)
    list=list(sec.order.par=rho,k0=k0star,threshold=u,tail.index=ti)
    list
  }
