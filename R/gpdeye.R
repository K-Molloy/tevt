#' Automated Approach for Interpreting the Hill-Plot
#'
#' An Implementation of the so called Eye-balling Technique proposed in Danielsson et al. (2016)
#' @param data vector of sample data
#' @param ws size of the moving window. \code{Default} is one percent of the data
#' @param epsilon size of the range in which the estimates can vary
#' @param h percentage of data inside the moving window that should lie in the tolerable range
#' @details The procedure searches for a stable region in the Hill-Plot by defining a moving window. Inside this window the estimates of the Hill estimator with respect to \code{k} have to be in a pre-defined range around the first estimate within this window. It is sufficient to claim that only \code{h} percent of the estimates within this window lie in this range. The smallest \code{k} that accomplishes this is then the optimal number of upper order statistics, i.e. data in the tail.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index by plugging in \code{k0} into the hill estimator}
#' @references Danielsson, J. and Ergun, L.M. and de Haan, L. and de Vries, C.G. (2016). Tail Index Estimation: Quantile Driven Threshold Selection.
#' @examples
#' data(danish)
#' gpd.eye(danish)
#' @export
gpd.eye = function(data,ws=0.01,epsilon=0.3,h=0.9){

    n=length(data)
    w=floor(ws*n)

    i=1:(n-1)
    x=sort(data,decreasing=TRUE)
    gamma=(cumsum(log(x[i]))/i)-log(x[i+1])
    alpha=1/gamma

    count=0;erg=c()
    for (k in 2:(length(alpha)-w)){
      for (i in 1:length(w)){
        if (alpha[k+i]<(alpha[k]+epsilon) && alpha[k+i]>(alpha[k]-epsilon)) {
          count=count+1
        } else {
          count=count
        }
      }
      erg[k]=count/w
    }
    erg=erg>h
    k0=min(which(erg==1))
    u=x[k0]
    ti=alpha[k0]
    list=list(k0=k0,threshold=u,tail.index=ti)
    list
  }
