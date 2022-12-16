##--------------------------------------------------##
## Plot variogram: adaptation of the function GeoModels::plot.GeoVariogram 

plot.GeoVariogram2 <- function (x, ...) 
{
  if (!inherits(x, "GeoVariogram")) 
    stop("Enter an object obtained from the function GeoVariogram\n")
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  ispatim = bivariate = FALSE
  if (!is.null(x$bint)) 
    ispatim = TRUE
  if (x$bivariate) 
    bivariate = TRUE
  lags = c(0, x$centers)
  numlags = length(lags)
  if (ispatim) 
    lagt = c(0, x$bint)
  else lagt = 0
  numlagt = length(lagt)
  slow = 0
  lags_m = seq(slow, max(x$centers), length.out = 150)
  if (ispatim) 
    lagt_m = seq(slow, max(x$bint), length.out = 150)
  else lagt_m = 0
  numlags_m = length(lags_m)
  numlagt_m = length(lagt_m)
  vario.main = "Spatial semi-variogram"
  vario.ylab = "Semi-Variogram"
  if (ispatim) {
    vario.main = "Space-time semi-variogram"
    vario.zlab = "Semi-Variogram"
  }
  if (bivariate) {
    par(mfrow = c(2, 2))
    plot.default(x$centers, x$variograms[1, ], main = "First semi-variogram", 
                 ylim = c(0, max(x$variograms[1, ])), xlim = c(0, 
                                                               max(x$centers)), xlab = "Distance", ylab = "Semi-Variogram", 
                 ...)
    if (min(x$variogramst) > 0) {
      ll1 = 0
      ll = max(x$variogramst)
    }
    if (min(x$variogramst) < 0) {
      ll1 = min(x$variogramst)
      ll = -min(x$variogramst)
    }
    plot.default(x$centers, x$variogramst, main = "Cross semi-variogram", 
                 ylim = c(ll1, ll), xlim = c(0, max(x$centers)), xlab = "Distance", 
                 ylab = "Semi-Variogram", ...)
    plot.default(x$centers, x$variogramst, main = "Cross semivariogram", 
                 ylim = c(ll1, ll), xlim = c(0, max(x$centers)), xlab = "Distance", 
                 ylab = "Semi-Variogram", ...)
    plot.default(x$centers, x$variograms[2, ], main = "Second semi-variogram", 
                 ylim = c(0, max(x$variograms[2, ])), xlim = c(0, 
                                                               max(x$centers)), xlab = "Distance", ylab = "Semi-Variogram", 
                 ...)
  }
  if (ispatim) {
    par(mfrow = c(2, 2), mai = c(0.7, 0.7, 0.3, 0.3), mgp = c(1.4,0.5, 0))
    plot.default(x$centers, x$variograms, xlab = expression(h), 
                 ylab = expression(gamma(h)), ylim = c(0, max(x$variograms)), 
                 xlim = c(0, max(x$centers)), main = "Marginal spatial semi-variogram")
    plot.default(x$bint, x$variogramt, xlab = expression(t), 
                 ylab = expression(gamma(t)), ylim = c(0, max(x$variogramt)), 
                 xlim = c(0, max(x$bint)), main = "Marginal temporal semi-variogram")
    evario = matrix(x$variogramst, nrow = length(x$centers), 
                    ncol = length(x$bint), byrow = TRUE)
    evario = rbind(c(0, x$variogramt), cbind(x$variograms, 
                                             evario))
    evario.grid = as.matrix(expand.grid(c(0, x$centers), 
                                        c(0, x$bint)))
    scatterplot3d::scatterplot3d(evario.grid[, 1], evario.grid[, 
                                                               2], c(evario), type = "h", highlight.3d = TRUE, 
                                 cex.axis = 0.7, cex.lab = 0.7, main = paste("Empirical", 
                                                                             vario.main), xlab = "Distance", ylab = "Time", 
                                 zlab = vario.zlab, mar = c(4, 3, 2, 2), mgp = c(0, 0, 0))
    par(mai = c(0.2, 0.3, 0.2, 0.2), mgp = c(1, 0.3, 0))
    persp(c(0, x$centers), c(0, x$bint), evario, xlab = expression(h), 
          ylab = expression(u), zlab = expression(gamma(h, u)), 
          ltheta = 90, shade = 0.75, ticktype = "detailed", 
          phi = 30, theta = 30, main = "Space-time semi-variogram", 
          cex.axis = 0.8, cex.lab = 0.8)
  }
  if (!ispatim && !bivariate) 
    plot.default(x$centers, x$variograms, ...)
  return(invisible())
}


