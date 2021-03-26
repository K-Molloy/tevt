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
