library(ppboot)
library(misc)
misc::clean_up()

mu   = .75
a    = .5
b    = 1.1
Tmax = 20
dt   = 1e-1

# Test error difference of fb and rb (they should be equal...)
set.seed(1234)
hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
N    = length(hp)

x = seq(0,Tmax,dt)
l = data.frame(x=x,y=l.hawkes(x,hp,mu,a,b))
L = data.frame(x=x,y=L.hawkes(x,hp,mu,a,b))

s  = interp.L(L,x=hp)

kern <- function(x,y){
  out = L.hawkes(x,y,mu,a,b)
}
# fb = pp.bootstrap1(L,s, NULL, FALSE, TRUE)
# rb = pp.bootstrap1(L,s, kern, FALSE, TRUE)


length(hp)

B = 399

Smax = max(L$y)



# system.time({
#   rb.boot1 = pp.boot(L,Smax,kern,B,NULL,TRUE)  
# })

system.time({
  rb.boot = pp.boot.par(L, Smax, kern, NULL, B, NULL)
})


# system.time({
#   fb.boot1 = pp.boot(L,Smax,NULL,B,TRUE)  
# })

system.time({
  fb.boot = pp.boot.par(L, Smax, NULL, B, NULL)
})

# parallel::stopCluster()


col1 = add.alpha('black',.95)
col2 = add.alpha('red',.75)
col3 = add.alpha('dodgerblue',.75)

# Check that fb and rb produce the same process with the same input Exp(1) variables and integrated conditional intensity
# layout(1)
# par(mar=c(3,3,1,1), bty='n', oma=c(0,0,0,0))
# plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n')
# 
# lines(L$x, L$y, col=col1)
# lines(l$x, l$y, col=col1)
# 
# points(hp, rep(0,N), pch=1, col=col1)
# points(rep(0,N), s, pch=16, col=col1)
# segments(hp,0,hp,s, lty=3, col=col1)
# segments(0,s,hp,s, lty=3, col=col1)
# 
# 
# points(fb, rep(0,N), pch=16, col=col2, cex=.75)
# segments(fb,0,fb,s, lty=3, col=col2)
# segments(0,s,fb,s, lty=3, col=col2)
# 
# points(rb, rep(0,N), pch=16, col=col3, cex=.75)
# segments(rb,0,rb,s, lty=3, col=col3)
# segments(0,s,rb,s, lty=3, col=col3)

hps = list()
for(i in 1:B){
  hp.tmp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  hps[[i]] = hp.tmp
}

layout(1)
# par(mfrow=c(1,3))

layout(matrix(c(1,2,3,4,4,4), nr=2, byrow=TRUE))

plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n', main="Hawkes")
for(i in 1:B){
  b.tmp = hps[[i]]
  L.tmp = L.hawkes(x,b.tmp, mu, a, b)
  lines(x,L.tmp, col=add.alpha(col1,.3))
}
lines(L$x, L$y, col=col1, lwd=2)
lines(L$x, L$y, col='white', lwd=2, lty=3)

plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n', main="FIB")
for(i in 1:B){
  b.tmp = fb.boot[[i]]
  L.tmp = L.hawkes(x,b.tmp, mu, a, b)
  lines(x,L.tmp, col=add.alpha(col2,.3))
}
lines(L$x, L$y, col=col1, lwd=2)

plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n', main="RB")
for(i in 1:B){
  b.tmp = rb.boot[[i]]
  L.tmp = L.hawkes(x,b.tmp, mu, a, b)
  lines(x,L.tmp, col=add.alpha(col3,.3))
}
lines(L$x, L$y, col=col1, lwd=2)

# unlist(lapply(fb.boot1,length))
# unlist(lapply(rb.boot1,length))

hist(unlist(lapply(hps,length)), border=NA, col=add.alpha(col1,.5), breaks=seq(0,100,1), prob=TRUE, ylim=c(0,.1), xlim=c(0,50))
par(new=TRUE)
hist(unlist(lapply(fb.boot,length)), border=NA, col=add.alpha(col2,.5), breaks=seq(0,100,1), prob=TRUE, ylim=c(0,.1), xlim=c(0,50))
par(new=TRUE)
hist(unlist(lapply(rb.boot,length)), border=NA, col=add.alpha(col3,.5), breaks=seq(0,100,1), prob=TRUE, ylim=c(0,.1), xlim=c(0,50))

quantile(unlist(lapply(hps,length)), probs = seq(0.75,.90,.01))
quantile(unlist(lapply(hps,length)))

