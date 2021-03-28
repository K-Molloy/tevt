#' A Double Bootstrap Procedure for Choosing the Optimal Sample Fraction
#'
#' An Implementation of the procedure proposed in Gomes et al. (2012) and Caeiro et al. (2016) for selecting the optimal sample fraction in tail index estimation.
#' @param data vector of sample data
#' @param B number of Bootstrap replications
#' @param epsilon gives the amount of the first resampling size \code{n1} by choosing \code{n1 = n^epsilon}. Default is set to \code{epsilon=0.995}
#' @details The Double Bootstrap procedure simulates the AMSE criterion of the Hill estimator using an auxiliary statistic. Minimizing this statistic gives a consistent estimator of the sample fraction \code{k/n} with \code{k} the optimal number of upper order statistics. This number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{second.order.par }{gives an estimation of the second order parameter \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail}
#' @references Gomes, M.I. and Figueiredo, F. and Neves, M.M. (2012). Adaptive estimation of heavy right tails: resampling-based methods in action. \emph{Extremes}, \bold{15}, 463--489.
#' @references Caeiro, F. and Gomes, I. (2016). Threshold selection in extreme value analysis. \emph{Extreme Value Modeling and Risk Analysis: Methods and Applications}, 69--86.
#' @examples
#' data(danish)
#' gpd.gomes(danish)
#' @export
gpd.gomes = function(data,B=1000,epsilon=0.995){
    n=length(data)
    n1=floor(n^epsilon)

    helphill=function (k,j) {
      xstat = sort(data, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^j)
      xihat
    }

    DAMSE=function(data){
      n=length(data)
      k=c(floor(n^0.995),floor(n^0.999))

      M11=helphill(k[1],1)
      M12=helphill(k[1],2)
      M13=helphill(k[1],3)

      M21=helphill(k[2],1)
      M22=helphill(k[2],2)
      M23=helphill(k[2],3)

      #tau=1
      W.k1.t1=(M11-sqrt(M12/2))/(sqrt(M12/2)-(M13/6)^(1/3))
      W.k2.t1=(M21-sqrt(M22/2))/(sqrt(M22/2)-(M23/6)^(1/3))
      rho.k1.t1=-abs(3*(W.k1.t1-1)/(W.k1.t1-3))
      rho.k2.t1=-abs(3*(W.k2.t1-1)/(W.k2.t1-3))

      #tau=0
      W.k1.t0=(log(M11)-0.5*log(M12/2))/(0.5*log(M12/2)-log(M13/6)/3)
      W.k2.t0=(log(M21)-0.5*log(M22/2))/(0.5*log(M22/2)-log(M23/6)/3)
      rho.k1.t0=-abs(3*(W.k1.t0-1)/(W.k1.t0-3))
      rho.k2.t0=-abs(3*(W.k2.t0-1)/(W.k2.t0-3))

      chi.t1=median(c(rho.k1.t1,rho.k2.t1))
      chi.t0=median(c(rho.k1.t0,rho.k2.t0))
      I.t1=c((rho.k1.t1-chi.t1)^2,(rho.k2.t1-chi.t1)^2)
      I.t0=c((rho.k1.t0-chi.t0)^2,(rho.k2.t0-chi.t0)^2)

      if (sum(I.t0)<=sum(I.t1)) {
        tau=0; rho=rho.k2.t0
      } else {
        tau=1; rho=rho.k2.t1
      }

      U=c();
      xstat = sort(data, decreasing = TRUE)
      for (i in 1:k[2]) {
        U[i]=i*(log(xstat[i]) - log(xstat[i + 1]))
      }

      i=1:k[2]
      dk=mean((i/k[2])^(-rho))

      Dk=function(a){
        D=mean((i/k[2])^(-a)*U[i])
      }

      beta=(k[2]/n)^rho*(dk*Dk(0)-Dk(rho))/(dk*Dk(rho)-Dk(2*rho))
      erg=c(beta,rho)
      erg
    }
    secorder=DAMSE(data)

    n2=floor((n1^2/n))+1
    khelp1=2:(n1-1)
    khelp2=2:(n2-1)
    tn1=matrix(nrow=B,ncol=n1-2)
    tn2=matrix(nrow=B,ncol=n2-2)
    for (l in 1:B){
      x2=sample(data,n2,replace=TRUE)
      x1=c(x2,sample(data,(n1-n2),replace=TRUE))
      x2=sort(x2,decreasing=TRUE)
      x1=sort(x1,decreasing=TRUE)
      i1=2:(n1-1)
      i2=2:(n2-1)
      m1=floor(i1/2)
      m2=floor(i2/2)
      Tkn1.1=(cumsum(log(x1[i1]))/i1)-log(x1[i1+1])
      Tkn1.2=(cumsum(log(x1[m1]))/i1)-log(x1[m1+1])
      Tkn2.1=(cumsum(log(x2[i2]))/i2)-log(x2[i2+1])
      Tkn2.2=(cumsum(log(x2[m2]))/i2)-log(x2[m2+1])
      tn1[l,]=Tkn1.2-Tkn1.1
      tn2[l,]=Tkn2.2-Tkn2.1
    }
    tn1=tn1^2
    tn2=tn2^2
    tn1star=colMeans(tn1)
    tn2star=colMeans(tn2)
    k1star=which.min(tn1star)
    k2star=which.min(tn2star)

    rho=secorder[2]
    Exp=2/(1-2*rho)
    k0star=floor((1-2^rho)^(Exp)*((k1star^2)/k2star))+1
    u=sort(data,decreasing=TRUE)[k0star]
    ti=1/helphill(k0star,1)
    list=list(sec.order.par=secorder,k0=k0star,threshold=u,tail.index=ti)
    list
  }
