#' QQ-Plot against the generalized Pareto distribution for given number of exceedances
#'
#' Plots the empirical observations above a given threshold against the theoretical quantiles of a generalized Pareto distribution.
#' @param data vector of sample data
#' @param nextremes number of exceedances
#' @param scale scale parameter of GPD
#' @param shape shape parameter of GPD
#' @details If the fitted GPD model provides a reasonable approximation of the underlying sample data the empirical and theoretical quantiles should coincide. So plotting them against each other should result in a straight line. Deviations from that line speak for a bad model fit and against a GPD assumption.
#' @return The straight red line gives the line of agreement. The dashed lines are simulated 95 percent confidence intervals. Therefor the fitted GPD model is simulated 1000 times using Monte Carlo. The sample size of each simulation equals the number of exceedances.
#' @examples
#' data=rexp(1000) # GPD with scale=1, shape=0
#' gpd.qq(data,1000,1,0) # Generate QQ plot
#' @export
gpd.qq =  function(data,nextremes,scale,shape){

    xstat=sort(data,decreasing=T)[1:nextremes]

    u=sort(data,decreasing=T)[nextremes]

    i=ppoints(nextremes)

    x=qgpd(i,loc=u,scale=scale,shape=shape)

    sim.new.quant=function(l){
      simdata=rgpd(nextremes,loc=u,scale=scale,shape=shape)
      sort(simdata)
    }

    sim.q=sapply(1:1000,sim.new.quant)

    ci.q=apply(sim.q,quantile,MARGIN=1,probs=c(0.025,0.975))

    make_qq_plot(x, u, xstat, ci.q)

}


make_qq_plot = function(x, u, xstat, ci.q){

  axislims=c(u,max(xstat))
  plot(x,sort(xstat),pch="x",xlab="theoretical quantiles",ylab="empirical quantiles",
       xlim=axislims,ylim=axislims)
  abline(0,1,col="red")
  lines(x,ci.q[1,],lty=2)
  lines(x,ci.q[2,],lty=2)
}

#' QQ-Estimator-Plot
#'
#' Plots the QQ-Estimator against the upper order statistics
#' @param data vector of sample data
#' @param kmin gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}
#' @param conf.int \code{logical}. If FALSE (default) no confidence intervals are plotted
#' @details The QQ-Estimator is a Tail Index Estimator based on regression diagnostics. Assuming a Pareto tail behaviour of the data at hand a QQ-Plot of the theoretical quantiles of an exponential distribution against the empirical quantiles of the log-data should lead to a straight line above some unknown upper order statistic \code{k}. The slope of this line is an estimator for the tail index. Computing this estimator via linear regression for every \code{k} the plot should stabilize for the correct number of upper order statistics, denoted \code{k0} here.
#' @return The plot shows the values of the QQ-Estimator with respect to \code{k}. See references for more information.
#' @references Kratz, M. and Resnick, S.I. (1996). The QQ-estimator and heavy tails. \emph{Stochastic Models}, \bold{12(4)}, 699--724.
#' @examples
#' data(danish)
#' gpd.qq.est(danish)
#' @export
gpd.qq.est = function(data,kmin=5,conf.int=FALSE){
    data=sort(data)
    n=length(data)

    QQplot=function(k){
      i=1:k
      xaxis=-log(1-(i/(k+1)))
      yaxis=log(data[n-k+i])
      SL=as.vector(lm(yaxis~xaxis)[[1]][2])
      return(1/SL)
    }
    x=1:n
    y=sapply(x,QQplot)
    plot(x[kmin:(n-1)],y[kmin:(n-1)],type="l",xlim=c(15,n),xlab="Order Statistics",ylab="qq-alpha")
    if (conf.int == TRUE){
      confint.d = 1/y - (qnorm(0.975) * sqrt(2)/(y*sqrt(x)))
      confint.u = 1/y + (qnorm(0.975) * sqrt(2)/(y*sqrt(x)))
      lines(x[kmin:(n-1)],1/confint.d[kmin:(n-1)], lty = 2, col = "blue")
      lines(x[kmin:(n-1)],1/confint.u[kmin:(n-1)], lty = 2, col = "blue")
    }
  }
