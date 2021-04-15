#' Handler function for a series of Bias-based procedures for Choosing the Optimal Sample Fraction
#'
#' @param data vector of sample data
#' @param method selected method \code{default} is set to \code{DK}
#' @param ... extra parameters for respective functions
#' @details Modifying the method parameter allows access to 4 similar bias based procedures
#' @return
#' \item{second.order.par }{gives an estimation of the second order parameter \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail}
#' @examples
#' data(danish)
#' gpd.OSF(danish)
#' @export
gpd.OSF = function(data, method="DK", ...){

  if (method=="DK"){
    DK(data, ...)
  } else if(method=="GH"){
    GH(data, ...)
  } else if(method=="RT"){
    RT(data, ...)
  } else if(method=="PS"){
    PS(data, ...)
  }
}
#' A Bias-based procedure for Choosing the Optimal Sample Fraction
#'
#' An Implementation of the procedure proposed in Drees & Kaufmann (1998) for selecting the optimal sample fraction in tail index estimation.
#' @param data vector of sample data
#' @param r tuning parameter for the stopping criterion. \code{default} is set to \code{1}. Change only if recommended by the output.
#' @details The procedure proposed in Drees & Kaufmann (1998) is based on bias reduction. A stopping criterion with respect to \code{k} is implemented to find the optimal tail fraction, i.e. \code{k/n} with \code{k} the optimal number of upper order statistics. This number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. If the above mentioned stopping criterion exceedes a certain value \code{r}, the bias of the assumed extreme model has become prominent and therefore \code{k} should not be chosen higher. For more information see references.
#' @return
#' \item{second.order.par }{gives an estimation of the second order parameter \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail}
#' @references Drees, H. and Kaufmann, E. (1998). Selecting the optimal sample fraction in univariate extreme value estimation. \emph{Stochastic Processes and their Applications}, \bold{75(2)}, 149--172.
#' @examples
#' data(danish)
#' gpd.OSF(danish, method="DK")
DK = function(data,r=1){
    n=length(data)

    gamma.tilde=helphill(data, 2*sqrt(n))
    r.n=r*2.5*gamma.tilde*n^0.25

    kn=function(k){
      i=1:k
      x=sort(data,decreasing=TRUE)
      h=(cumsum(log(x[i]))/i)-log(x[i+1])
      h2=sqrt(i)*abs(h-h[k])
      Max=max(which(h2==max(h2)))
      Max
    }
    khelp=sapply(1:(n-1),kn)
    khelp2=khelp>r.n
    kbar_r.n=min(which(khelp2==1))

    if (kbar_r.n==Inf) {
      print("Warning: no k_n>r_n found. Use smaller r_n. Try r=0.9")
    } else {
      xi=0.7
      lambda=0.6
      r.n.xi=r.n^xi

      khelp3=khelp>r.n.xi
      kbar_r.n.xi=min(which(khelp3==1))

      Z=log(kn(floor(lambda*kbar_r.n.xi)))-log(kn(kbar_r.n.xi))
      rho=-(Z/log(lambda))+0.5

      if (rho>=0) {
        print("Warning: rho>0. Method fails here.")
      } else {
        Exp=1/(1-2*rho)
        kk=(kbar_r.n.xi/(kbar_r.n^xi))^(1/(1-xi))
        k0star=(((1-2*rho)^(1/rho))*((-2*rho*gamma.tilde^2)^Exp)*kk)
        k0star=floor(k0star)
        ti=1/helphill(data, k0star)
        u=sort(data,decreasing=TRUE)[k0star]
        list=list(sec.order.par=rho,k0=k0star,threshold=u,tail.index=ti)
        list
      }
    }
}

#' A Bias-based procedure for Choosing the Optimal Threshold
#'
#' An Implementation of the procedure proposed in Guillou & Hall(2001) for selecting the optimal threshold in extreme value analysis.
#' @param data vector of sample data
#' @details The procedure proposed in Guillou & Hall (2001) is based on bias reduction. Due to the fact that the log-spacings of the order statistics are approximately exponentially distributed if the tail of the underlying distribution follows a Pareto distribution, an auxilliary statistic with respect to \code{k} is implemented with the same properties. The method then behaves like an asymptotic test for mean \code{0}. If some critical value \code{crit} is exceeded the hypothesis of zero mean is rejected. Thus the bias has become too large and the assumed exponentiality and therefore the assumed Pareto tail can not be hold.  From this an optimal number of \code{k} can be found such that the critical value is not exceeded. This optimal number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by
#' coosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Guillou, A. and Hall, P. (2001). A Diagnostic for Selecting the Threshold in Extreme Value Analysis. \emph{Journal of the Royal Statistical Society}, \bold{63(2)}, 293--305.
#' @examples
#' data(danish)
#' gpd.OSF(danish, method="GH")
GH = function(data){
    n=length(data)

    Qn=function(k){

      Tn=function(k){
        U=c();help=c();
        xstat = sort(data, decreasing = TRUE)
        for (i in 1:k) {
          U[i]=i*(log(xstat[i]) - log(xstat[i + 1]))
          help[i]=(k-2*i+1)*U[i]
        }
        Z=sum(help)
        N=mean(U)
        T=sqrt(3/k^3)*(Z/N)
        T
      }

      start=k-floor(k/2)
      end=k+floor(k/2)
      y=start:end
      x=sapply(y,Tn)
      erg=1/(2*floor(k/2)+1)*sum(x^2)
      erg2=sqrt(erg)
      crit=1.25
      as.numeric(erg2>=crit)
    }

    kmax=floor(n/1.5)
    i=1
    while (Qn(i)==0){
      i=i+1
      if (i==kmax) break
    }

    u=sort(data,decreasing=TRUE)[i]

    ti=1/helphill(data, i)
    list=list(k0=i,threshold=u,tail.index=ti)
    list
}

#' Sample Path Stability Algorithm
#'
#' An Implementation of the heuristic algorithm for choosing the optimal sample fraction proposed in Caeiro & Gomes (2016), among others.
#' @param data vector of sample data
#' @param j digits to round to. Should be \code{0} or \code{1} (default)
#' @details The algorithm searches for a stable region of the sample path, i.e. the plot of a tail index estimator with respect to \code{k}. This is done in two steps. First the estimation of the tail index for every \code{k} is rounded to \code{j} digits and the longest set of equal consecutive values is chosen. For this set the estimates are rounded to \cite{j+2} digits and the mode of this subset is determined. The corresponding biggest k-value, denoted \code{k0} here, is the optimal number of data in the tail.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Caeiro, J. and Gomes, M.I. (2016). Threshold selection in extreme value analysis. \emph{Extreme Value Modeling and Risk Analysis:Methids and Applications}, 69--86.
#' @references Gomes, M.I. and Henriques-Rodrigues, L. and Fraga Alves, M.I. and Manjunath, B. (2013). Adaptive PORT-MVRB estimation: an empirical comparison of two heuristic algorithms. \emph{Journal of Statistical Computation and Simulation}, \bold{83}, 1129--1144.
#' @references Gomes, M.I. and Henriques-Rodrigues, L. and Miranda, M.C. (2011). Reduced-bias location-invariant extreme value index estimation: a simulation study. \emph{Communications in Statistic-Simulation and Computation}, \bold{40}, 424--447.
#' @examples
#' data(danish)
#' gpd.OSF(danish, method="PS")
PS = function(data,j=1){
    n=length(data)
    x1=sort(data, decreasing=TRUE)
    i=1:(n-1)
    h=(cumsum(log(x1[i]))/i)-log(x1[i+1])

    help=round(h,j)
    l=rle(help)$lengths
    Max=which(l==max(l))
    kmin=sum(l[1:(Max-1)])+1
    kmax=kmin+max(l)-1
    run=kmax-kmin

    K=round(h[kmin:kmax],j+2)

    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    if (sum(duplicated(K))==0) {
      modus=K[length(K)]
    } else{
      modus=Mode(K)
    }

    k0star=max(which(round(h,j+2)==modus))
    u=x1[k0star]
    ti=1/h[k0star]
    list=list(k0=k0star,threshold=u,tail.index=ti)
    list
}

#' Adaptive choice of the optimal sample fraction in tail index estimation
#'
#' An implementation of the minimization criterion proposed in Reiss & Thomas (2007).
#' @param data vector of sample data
#' @param beta a factor for weighting the expression below. Default is set to \code{beta=0}
#' @param kmin gives a minimum value for \code{k}. Default ist set to \code{kmin=2}
#' @details The procedure proposed in Reiss & Thomas (2007) chooses the lowest upper order statistic \code{k} to minimize the expression
#' \code{1/k sum_i=1^k i^beta |gamma_i-median(gamma_1,...,gamma_k)|}
#' or an alternative of that by replacing the absolute deviation with a squared deviation and the median just with \code{gamma_k}, where \code{gamma} denotes the Hill estimator
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail for both metrics, i.e. the absolute and squared deviation.}
#' \item{threshold}{the corresponding thresholds.}
#' \item{tail.index}{the corresponding tail indices}
#' @references Reiss, R.-D. and Thomas, M. (2007). Statistical Analysis of Extreme Values: With Applications to Insurance, Finance, Hydrology and Other Fields. \emph{Birkhauser, Boston}.
#' @examples
#' data(danish)
#' gpd.OSF(danish, method="RT")
RT = function(data,beta=0,kmin=2){
    n=length(data)

    x1=sort(data, decreasing=TRUE)
    i=1:(n-1)
    j=i^beta
    h=(cumsum(log(x1[i]))/i)-log(x1[i+1])

    med=c();
    for (i in 1:(n-1)){
      med[i]=median(h[1:i])
    }
    erg1=c(); erg2=c()
    for (k in 1:(n-1)){
      erg1[k]=sum(j[1:k]*abs(h[1:k]-med[k]))/k
      erg2[k]=sum(j[1:k]*(h[1:k]-h[k])^2)/k
    }
    rt1=which.min(erg1[kmin:(n-1)])
    rt2=which.min(erg2[kmin:(n-1)])
    list=list(k0=c(rt1,rt2),threshold=c(x1[rt1],x1[rt2]),tail.index=c(1/h[rt1],1/h[rt2]))
    list
  }

