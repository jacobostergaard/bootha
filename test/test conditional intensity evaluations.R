library(ppboot)
library(misc)
misc::clean_up()

set.seed(1234)

Tmax = 1000
mu   = .75
a    = .5
b    = 1.1
dt   = 1e-2

hp = unlist(simulateHawkes(mu,a,b,Tmax))


# Check that evaluations work when x is not a full obsevation interval
  y = hp[1:10]

  x1 = seq(0,10,.1)
  x2 = seq(5,10,.1)
  x3 = seq(0,5,.1)
  x4 = seq(2.5,7.4,.1)
  
  par(mfrow = c(2,1), mar=c(3,3,1,1), oma=c(1,1,1,1), bty='n')
  y1 = l.hawkes(x1,y)
  y2 = l.hawkes(x2,y)
  y3 = l.hawkes(x3,y)
  y4 = l.hawkes(x4,y)
  
  plot(x1,y1, type='l', lwd=5)
  lines(x2,y2, type='l', col=add.alpha('red',.8), lwd=3)
  lines(x3,y3, type='l', col=add.alpha('blue',.5), lwd=3)
  lines(x4,y4, type='l', col=add.alpha('orange',.5), lwd=3)
  
  y1 = L.hawkes(x1,y)
  y2 = L.hawkes(x2,y)
  y3 = L.hawkes(x3,y)
  y4 = L.hawkes(x4,y)
  
  plot(x1,y1, type='l', lwd=5)
  lines(x2,y2, type='l', col=add.alpha('red',.8), lwd=3)
  lines(x3,y3, type='l', col=add.alpha('dodgerblue',.5), lwd=3)
  lines(x4,y4, type='l', col=add.alpha('orange',.5), lwd=3)
  





# Check that transformations correspond to correct evaluations
  x = seq(0,1000,.1)
  y = l.hawkes(x,hp)
  Y = L.hawkes(x,hp)
  
  s = interp.L(cbind(x,Y),hp)
  n = length(hp)
  
  par(mfrow = c(2,1), mar=c(3,3,1,1), oma=c(1,1,1,1), bty='n')
  plot(x,y, type='l', lwd=1, xlim=c(25,50), ylim=c(0,4))
  points(hp,rep(0,n))
  
  plot(x,Y, type='l', lwd=1, xlim=c(25,50),ylim=c(40,75))
  points(hp,rep(40,n))
  points(rep(25,n),s)
  segments(25,s,hp,s, lty=3)
  segments(hp,40,hp,s, lty=3)

  # MSE when inverting L from s-time to t-time:
  sum((hp-interp.L(cbind(x,Y),y=s))^2)
  

# Check that transformed waiting times are Exp(1)
  layout(1)
  w = diff(s)
  hist(w, breaks=seq(0,1000,.1), xlim=c(0,10), prob=TRUE)
  curve(dexp,add=TRUE, col='red')

  
  
  
