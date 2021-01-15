library(misc)
library(bootha)

set.seed(1234)
Tmax = 15
pars = c(1, 1.5/2, 2)
s = simulateHawkes(Tmax, pars, verbose=TRUE)
x = seq(0,Tmax,length=5000)

misc::set_cobalt_plot_bg()
rmv = (s$x %in% setdiff(s$x,s$y))
plot(s$x[rmv], s$u[rmv], col=add.alpha('red',.9), pch=16, xlim=c(0,Tmax), ylim=c(0,max(pretty(s$m))))
segments(s$x[!rmv],0,s$x[!rmv],s$u[!rmv], lty=3, lwd=1, col=add.alpha('white',.5))
points(s$x[!rmv], s$u[!rmv], col=add.alpha('dodgerblue',.8), pch=16)
points(s$y, rep(0,length(s$y)), pch=124, col='white')
lines(x,intensity(x,s$y,pars), type='s', col=add.alpha('white',.5))
# lines(c(0,s$x),c(s$m,s$m[length(s$m)]), type='s', col='orange')
segments(c(0,s$x[-length(s$x)]),s$m,s$x,s$m, lwd=2, lty=1, col=add.alpha('orange',.8))

