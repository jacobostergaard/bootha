library(bootha)
library(misc)
clean_up()
set.seed(1234)

pow_kernel <- function(x,a=.5,b=.75,k=.8){
  a*b/( (1+b*x)^(1+k) )
  
}
exp_kernel <- function(x,a=.5,b=.75){
  a*exp(-b*x)
}

a = .5
b = .75
k = .8

x = seq(0,10,.1)


layout(1)
plot(x,pow_kernel(x, a,b,k), type='l')
lines(x, exp_kernel(x,a,b), lty=3)

integrate(exp_kernel, lower = 0, upper = Inf)$value*b/a
integrate(pow_kernel, lower = 0, upper = Inf)$value*k/a


x = seq(0,10,.1)
k = 1.33
b = 10/250.05
a = 1*(b^k)

k = 0.5
a = .25
b = 1.1
B = 1.01
B = 1
m = 9e4
plot(x,m^B*exp_kernel(x,a,b), lty=3, type='l')
lines(x,m^B*pow_kernel(x,a,b,k))



Tmax = 3000
pars = c(.05,.8,.2,.9)
pars[2]/pars[4]
system.time({
  pp = bootha::simulateHawkes(Tmax, pars)  
})

x  = seq(0,Tmax,.1)
y  = intensity(x, pp, pars)
Nt = length(pp)


layout(1:2)
plot(x,y, type='l', ylim=c(-.1*max(y),max(y)))
points(pp, rep(-.1*max(y),Nt), pch=124,col = add.alpha('black',.5))
plot(x,integrated_intensity(x, pp, pars), type='l')

system.time({
  res = mle(pp, Tmax, pars, dt=1) # fit power kernel
  res2 = mle(pp, Tmax, pars[1:3], dt=1, ui0=diag(c(1,-1,1)), ci0=c(0,-1,0))   # fit exp kernel
})

Nt
pars
res
res2
res[2]/res[4]
pars[2]/pars[4]
