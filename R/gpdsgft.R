#' Sequential Goodness of Fit Testing for the Generalized Pareto Distribution
#'
#' An implementation of the sequential testing procedure proposed in Thompson et al. (2009) for automated threshold selection
#' @param data vector of sample data
#' @param thresholds a sequence of pre-defined thresholds to check for GPD assumption
#' @details The procedure proposed in Thompson et al. (2009) is based on sequential goodness of fit testing. First, one has to choose a equally spaced grid of posssible thresholds. The authors recommend 100 thresholds between the 50 percent and 98 percent quantile of the data, provided there are enough observations left (about 100 observations above the last pre-defined threshold). Then the parameters of a GPD for each threshold are estimated. One can show that the differences of subsequent scale parameters are approximately normal distributed. So a Pearson chi-squared test for normality is applied to all the differences, striking the smallest thresholds out until the test is not rejected anymore.
#' @return
#' \item{threshold}{the threshold used for the test}
#' \item{num.above}{the number of observations above the given threshold}
#' \item{p.values}{raw p-values for the thresholds tested}
#' \item{ForwardStop}{transformed p-values according to the ForwardStop criterion. See G'Sell et al (2016) for more information}
#' \item{StrongStop}{transformed p-values according to the StrongStop criterion. See G'Sell et al (2016) for more information}
#' \item{est.scale}{estimated scale parameter for the given threshold}
#' \item{est.shape}{estimated shape parameter for the given threshold}
#' @references Thompson, P. and Cai, Y. and Reeve, D. (2009). Automated threshold selection methods for extreme wave analysis. \emph{Coastal Engineering}, \bold{56}(10), 1013--1021.
#' @references G'Sell, M.G. and Wager, S. and Chouldechova, A. and Tibshirani, R. (2016). Sequential selection procedures and false discovery rate control. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} \bold{78}(2), 423--444.
#' @examples
#' data=rexp(1000)
#' u=seq(quantile(data,.1),quantile(data,.9),,100)
#' A=gpd.sgft(data,u);A
#' @export
gpd.sgft = function(data,thresholds){
    u=thresholds

    pearson.test=function (x, n.classes = ceiling(2 * (n^(2/5))), adjust = TRUE)
    {
      DNAME <- deparse(substitute(x))
      x <- x[complete.cases(x)]
      n <- length(x)
      if (adjust) {
        dfd <- 2
      }
      else {
        dfd <- 0
      }
      num <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
      count <- tabulate(num, n.classes)
      prob <- rep(1/n.classes, n.classes)
      xpec <- n * prob
      h <- ((count - xpec)^2)/xpec
      P <- sum(h)
      pvalue <- pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
      RVAL <- list(statistic = c(P = P), p.value = pvalue, method = "Pearson chi-square normality test",
                   data.name = DNAME, n.classes = n.classes, df = n.classes -
                     1 - dfd)
      class(RVAL) <- "htest"
      return(RVAL)
    }

    shape=c();scale=c();exc=c()
    for (l in 1:length(u)) {
      est=gpd.fit(data,threshold=u[l])
      exc[l]=est$n.exceed
      shape[l]=est[[4]][2,1]
      scale[l]=est[[4]][1,1]
      rm(est)
    }

    tu=c();
    for (j in 1:length(u)) {
      tu[j]=scale[j]-shape[j]*u[j]
    }
    t=diff(tu)

    p.v=c();
    for (q in 1:(length(t)-2)){
      p.v[q]=pearson.test(t[q:length(t)])$p.value
    }
    stops=pSeqStop(p.v)
    names=c("testnum","threshold","num.above","p.values",
            "ForwardStop","StrongStop","est.scale","est.shape")
    erg=matrix(ncol=8,nrow=length(u)-3)
    colnames(erg)=names
    erg[,1]=1:(length(u)-3)
    erg[,2]=u[1:(length(u)-3)]
    erg[,3]=exc[1:(length(u)-3)]
    erg[,4]=c(p.v)
    erg[,5]=c(stops$ForwardStop)
    erg[,6]=c(stops$StrongStop)
    erg[,7]=scale[1:(length(u)-3)]
    erg[,8]=shape[1:(length(u)-3)]

    erg
  }
