# Functions for the Fixed intensity bootstrap (FIB)

fib <- function(L, type='par', plotit=FALSE){
  # Fixed intensity bootstrap, currently only parametric bootstrap available
  # Here L is the integrated conditional intensity in (x,y) format
  
  # Algorithm:
  #   1. Sample waiting times in scaled time domain
  #   2. Invert to observed time domain using input integrated intensity
  
  Smax = max(L$y)
  n    = ceiling(Smax*2)      # approximate counts in [0,T] is lam*T, hence twice the expected counts, should avoid the while loop.
  W    = parboot_wt(n)   # waiting times
  
  while(sum(W)<Tmax){       # If we still haven't reached the end of the interval [0,T] yet, simulate more waiting times
    W = c(W,parboot_wt(1))
  }
  
  s = cumsum(W)   # event times
  s = s[s<Smax]   # truncate to interval [0,T]
  
  y.new = s
  x.new = interp.L(L,y = y.new)
  
  if(plotit){
    plot(L, type='l', lwd=1.5, bty='n')
    points(rep(0,length(y.new)),y.new, pch=16, cex=.75)
    points(x.new,rep(0,length(x.new)), pch=6, cex=.75)
    segments(0,y.new,x.new,y.new, lty=3)
    segments(x.new,0,x.new,y.new, lty=3)    
  }

  out = data.frame(t = x.new, s = y.new)
  return(out)
  
}

# fib(L, plotit=TRUE)

# layout(1)
# par(mar=c(3,3,1,1))
# fib(L, plotit=TRUE)
# 
# x = 1:1000
# y = rep(1,1000)
# L = data.frame(x=x, y=cumsum(y)/1e2)
# 
# x = 1:100
# y1 = 3*x^2
# y2 = .75*x^2
# y3 = .5*x^2
# y4 = 1.1*x^2
# y5 = rep(1e3,100)
# 
# y = c(y1,y2,y3,y4,y5)/1e5
# x = 1:length(y)
# L = data.frame(x=x, y=cumsum(y))
# layout(1:2)
# plot(x,y, type='l')
# plot(L, type='l')


