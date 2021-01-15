library(ppboot)
library(misc)
misc::clean_up()

mu   = .75
a    = .5
b    = 1.1
Tmax = 100
dt   = 1e-2

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

fb = pp.bootstrap1(L,s, NULL, FALSE, TRUE)
rb = pp.bootstrap1(L,s, kern, FALSE, TRUE)


B = 100

Smax = max(L$y)

system.time({
  rb.boot.nopar = pp.boot(L,Smax,kern,B,NULL,parallel = FALSE, verbose=TRUE)
})

system.time({
  rb.boot.par = pp.boot(L,Smax,kern,B,NULL,parallel = TRUE)
})


system.time({
  fb.boot.nopar = pp.boot(L,Smax,NULL,B,NULL,parallel = FALSE, verbose=TRUE)
})

system.time({
  fb.boot.par = pp.boot(L,Smax,NULL,B,NULL,parallel = TRUE)
})

col1 = add.alpha('black',.95)
col2 = add.alpha('red',.75)
col3 = add.alpha('dodgerblue',.75)

x = seq(0,Tmax,.1)
par(mfrow=c(2,2))

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200))
for(i in 1:B){
  y = L_exp_kernel(x, rb.boot.nopar[[i]], mu, a, b)
  lines(x,y, col=col2)
}

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200))
for(i in 1:B){
  y = L_exp_kernel(x, rb.boot.par[[i]], mu, a, b)
  lines(x,y, col=col2)
}

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200))
for(i in 1:B){
  y = L_exp_kernel(x, fb.boot.nopar[[i]], mu, a, b)
  lines(x,y, col=col3)
}

plot(0,0, type='n', xlim=c(0,100), ylim=c(0,200))
for(i in 1:B){
  y = L_exp_kernel(x, fb.boot.par[[i]], mu, a, b)
  lines(x,y, col=col3)
}

