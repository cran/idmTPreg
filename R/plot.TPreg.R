plot.TPreg <-
function(x, covar, rug=TRUE,  main, ylab, xlab, Ylim,  ...)
{
  if (missing(xlab)) xlab="time"
  if (x$transition=="all") {
    nfun <- 4
    nam.p <- c("11", "12", "13", "23")
    
  }   
  else {
    nfun <- 1
    nam.p <- x$transition
  }
  if (missing(covar)) covar <- colnames(x[[1]]$coefficients)[-1]
  if (!is.character(covar)) stop(" 'covar' must be a single character string referring the name of a covariable")
  ncovar <- length(covar)
  n.plots <- nfun*ncovar  
  mar0=c(4,3.9,4,3.9)
  if (n.plots == 0) 
    stop("No terms to plot - nothing for plot.TPreg() to do.")
  
  ppp <- n.plots
  c <- r <- trunc(sqrt(ppp))
  if (c < 1)       r <- c <- 1
  if (c * r < ppp) c <- c + 1
  if (c * r < ppp) r <- r + 1
  oldpar <- par(mfrow = c(r, c), mar = mar0)
  mis.lab <- missing(ylab) 
  mis.main <- missing(main)
  mis.Ylim <- missing(Ylim)
  for (i.fun in 1:nfun){
    for (i.covar in 1:ncovar){
      if (mis.Ylim) {ylim = NULL
      }
      else {ylim = Ylim[[i.covar]]
      }
      col <- covar[i.covar]
      time11 <- x[[i.fun]]$time
      i.co <- x[[i.fun]]
      co11 <- i.co$coeff
      if (mis.lab) ylab <- paste("effect", "of", col)
      if (mis.main) main <- paste("TP:", nam.p[i.fun] )
      plot(time11, co11[ , col], type = "s", xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
      lines(time11, i.co$LWL[,col], type = "s", col = "lightseagreen")
      lines(time11, i.co$UPL[,col], type = "s", col = "lightseagreen")
      lines(time11, rep(0,length(time11)), col = "red", lty = 2)
      if (rug) rug(as.numeric(time11 ), ...)
    }
  }
}
