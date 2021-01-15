add_alpha <- function(col, alpha=1){
    ## Add an alpha value to a color
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
            rgb(x[1], x[2], x[3], alpha=alpha))
}


make_ci <- function(x){
  out = c(mean(x)-1.96*sd(x),mean(x)+1.96*sd(x))
  return(out)
}


ks_plot <- function(z, a = .05, col1='tomato2', col2='steelblue', axlabs = TRUE){
  # Transform rescaled event times into U(0,1) variables and plot against identity
  # Input:
  #   z:      rescaled spike times
  #   a:      confidence level for bounds
  #   col1:   color for the U(0,1) variables
  #   col2:   color for the confidence region
  #   axlabs: if TRUE, then labels both axes
  
  z = sort(z)
  N = length(z)
  u = 1-exp(-z)
  r = (1:N-.5)/N
  
  # Get critical values of Kolmogorov-Smirnov test
  # see: https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test or
  # "The Art of Computer Programming, Volume 2" by DE Knuth, Eq (15) in section 3.3.1
  crit_val = sqrt(-.5*log(a/2))
  
  x = seq(-1,2,.1)
  
  plot(0,0, type='n', xlim=c(0,1), ylim=c(0,1), bty='n', ann=FALSE)
  polygon(c(x,rev(x)), c(x+crit_val/sqrt(N),rev(x-crit_val/sqrt(N))), col=col2, border=NA)
  
  abline(0,1,col=add.alpha('black', .60), lwd=1.5)
  lines(r,u, col=add.alpha(col1, .75), lwd=1.5)
  
  if(axlabs){
    mtext(text = "Theoretical", side=1, line= 2.5, cex=.75)
    mtext(text = "Empirical", side=2, line=2.5, las=0, cex=.75)
  }
  
}




plot_acf <- function(x, lag=25, col1 = 'black', col2 = 'gray60', ...){
  
  ci  = 1.96/sqrt(length(x))
  tmp = acf(x, lag.max = lag, plot=FALSE)
  
  # Omit first obs. as this is 1 and therefore more difficult to see significant correlations
  x = tmp$lag[-1]
  y = tmp$acf[-1]
  
  plot(x,y, ylim=c(-2,2)*ci, type='h', col=col1, ...)
  abline(h=c(-1,1)*ci,lty=3, col=col2)
  
}


bin <- function(pp, dt){
  minmax = range(pp)
  n   = diff(minmax)
  vec = numeric(n+1)
  vec[pp] = 1
  vec = c(vec, rep(NA,14-n%%dt-1))
  tmp = matrix(vec, nr = ceiling(n/dt), nc=dt, byrow=TRUE)
  out = apply(tmp, 1, function(x) sum(x, na.rm = TRUE))  
  return(out)  
}

