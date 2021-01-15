library(ppboot)
library(misc)
misc::clean_up()

mu   = .75
a    = .5
b    = 1.1
Tmax = 100
dt   = 1e-1

# Test error difference of fb and rb (they should be equal...)
set.seed(1234)
hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
N    = length(hp)

x = seq(0,Tmax,dt)

l = data.frame(x=x,y=l_exp_kernel(x,hp,mu,a,b))
L = data.frame(x=x,y=L_exp_kernel(x,hp,mu,a,b))


s  = interpol(as.matrix(L),x=hp)

kern <- function(x,y){
  if(is.null(y)){
    y = numeric(0)
  }
  out = L_exp_kernel(x,y,mu,a,b)
}


B    = 399
Smax = max(L$y)

system.time({
  rb.boot.nopar = pp.boot(L,Smax,kern,B,s,parallel = TRUE)
})

system.time({
  rb.boot.par = pp.boot(L,Smax,kern,B,NULL,parallel = TRUE)
})


system.time({
  fb.boot.nopar = pp.boot(L,Smax,NULL,B,s,parallel = TRUE)
})

system.time({
  fb.boot.par = pp.boot(L,Smax,NULL,B,NULL,parallel = TRUE)
})

col1 = add.alpha('black',.5)
col2 = add.alpha('red',.2)
col3 = add.alpha('dodgerblue',.2)



par(mfrow=c(2,2), mar=c(4,4,1,1)+.5, oma=c(0,0,0,0))

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200), main="RB parametric", xlab="t", ylab=expression(Lambda(t)))
abline(h=seq(0,200,25), lty=3, col=col1)
for(i in 1:B){
  y = L_exp_kernel(x, rb.boot.par[[i]], mu, a, b)
  lines(x,y, col=col2)
}


plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200), main="RB non-parametric", xlab="t", ylab=expression(Lambda(t)))
abline(h=seq(0,200,25), lty=3, col=col1)
for(i in 1:B){
  y = L_exp_kernel(x, rb.boot.nopar[[i]], mu, a, b)
  lines(x,y, col=col2)
}


plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200), main="FIB parametric", xlab="t", ylab=expression(Lambda(t)))
abline(h=seq(0,200,25), lty=3, col=col1)
for(i in 1:B){
  y = L_exp_kernel(x, fb.boot.par[[i]], mu, a, b)
  lines(x,y, col=col3)
}

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200), main="FIB non-parametric", xlab="t", ylab=expression(Lambda(t)))
abline(h=seq(0,200,25), lty=3, col=col1)
for(i in 1:B){
  y = L_exp_kernel(x, fb.boot.nopar[[i]], mu, a, b)
  lines(x,y, col=col3)
}

notify("Bootstrap done!")


