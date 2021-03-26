#' GPD Accessory Plots
#'
#' The plot method plot.gpd provides four different plots for assessing fitted GPD model.
#'
#' 1 ~ Excess Distribution
#' 2 ~ Tail of Underlying Distribution
#' 3 ~ Scatterplot of Residuals
#' 4 ~ QQplot of Residuals
#'
#'
#' @param x GPD Object
#' @param pick Plot selection type
#' @param optlog optional argument for plots 1 and 2 giving a particular choice of logarithmic axes: "x" x-axis only; "y" y-axis only; "xy" both axes; "" neither axis.
#' @param extend optional argument for plots 1 and 2 expressing how far x-axis should extend as a multiple of the largest data value. This argument must take values greater than 1 and is useful for showing estimated quantiles beyond data.
#' @param labels optional argument for plots 1 and 2 specifying whether or not axes should be labelled
#' @param ... other graphics parameters
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' data(danish)
#' out = gpd.fit(danish, 10)
#' gpd.plot(out, pick=1)

#' @export
gpd.plot = function (x, pick=1, optlog = NA, extend = 1.5, labels = TRUE, ...)
{
  data = as.numeric(x$exceedances)
  threshold = x$threshold
  xi = x$par.ests[2]
  beta = x$par.ests[1]

  lastcurve = NULL

  if (pick >= 3) {
    excess = data - threshold
    res = logb(1 + (xi * excess)/beta)/xi
    lastcurve = NULL
  }
  if (pick == 3) {
    plot(res, ylab = "Residuals", xlab = "Ordering",
         ...)
    lines(lowess(1:length(res), res))
  }
  if (pick == 4)
    qplot(res, ...)
  if (pick == 1 || pick == 2) {
    plotmin = threshold
    if (extend <= 1)
      stop("extend must be > 1")
    plotmax = max(data) * extend
    xx = seq(from = 0, to = 1, length = 1000)
    z = qgpd(xx, xi, threshold, beta)
    z = pmax(pmin(z, plotmax), plotmin)
    ypoints = ppoints(sort(data))
    y = pgpd(z, xi, threshold, beta)
  }
  if (pick == 1) {
    type = "eplot"
    if (!is.na(optlog))
      alog = optlog
    else alog = "x"
    if (alog == "xy")
      stop("Double log plot of Fu(x-u) does\nnot make much sense")
    yylab = "Fu(x-u)"
    shape = xi
    scale = beta
    location = threshold
  }
  if (pick == 2) {
    type = "tail"
    if (!is.na(optlog))
      alog = optlog
    else alog = "xy"
    prob = x$p.less.thresh
    ypoints = (1 - prob) * (1 - ypoints)
    y = (1 - prob) * (1 - y)
    yylab = "1-F(x)"
    shape = xi
    scale = beta * (1 - prob)^xi
    location = threshold - (scale * ((1 - prob)^(-xi) -
                                       1))/xi
  }
  if (pick == 1 || pick == 2) {
    plot(sort(data), ypoints, xlim = range(plotmin, plotmax),
         ylim = range(ypoints, y, na.rm = TRUE), xlab = "",
         ylab = "", log = alog, axes = TRUE)
    lines(z[y >= 0], y[y >= 0])
    if (labels) {
      xxlab = "x"
      if (alog == "x" || alog == "xy" ||
          alog == "yx")
        xxlab = paste(xxlab, "(on log scale)")
      if (alog == "xy" || alog == "yx" ||
          alog == "y")
        yylab = paste(yylab, "(on log scale)")
      title(xlab = xxlab, ylab = yylab)
    }
    details = paste("threshold = ", format(signif(threshold,
                                                  3)), "   xi = ", format(signif(shape, 3)),
                    "   scale = ", format(signif(scale, 3)),
                    "   location = ", format(signif(location,
                                                    3)), sep = "")
    print(details)
    lastcurve = list(lastfit = x, type = type, dist = "gpd",
                     plotmin = plotmin, plotmax = plotmax, alog = alog,
                     location = as.numeric(location), shape = as.numeric(shape),
                     scale = as.numeric(scale))
  }
  invisible(lastcurve)
}
