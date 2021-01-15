library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts5.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax

mle  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)
O   = mle

x = seq(0,Tmax,1)
L = integrated_intensity(x, pp, mle)
s = interpol(cbind(x,L), x=pp)

s.obs = s
n  = 400
v  = rexp(n)

s  = cumsum(v)
w  = numeric(n)
tt = cumsum(w)

f <- function(x,O,tk,Sk){
  O[1]*(x-tk)+O[2]*(1-exp(-O[3]*(x-tk) ) )*Sk
}
dfdx <- function(x,O,tk,Sk){
  O[1]+O[2]*Sk*exp(-O[3]*(x-tk))
}
i=1
# for(i in 1:n){
#   k  = max(which(s[i]>s.obs))
#   Sk = sum(exp(-O[3]*(pp[k]-pp[1:k])))
#   tk = pp[k]
#   sk = s.obs[k]
#   g  = function(x) f(x, O, tk, Sk)+sk-s[i]
#   tt[i] = uniroot(g, c(0, 1e6))$root
# }

for(i in 1:n){
  k  = max(which(s[i]>s.obs))
  Sk = sum(exp(-O[3]*(pp[k]-pp[1:k])))
  tk = pp[k]
  sk = s.obs[k]
  
  if(i>1){
    x0 = tt[i-1]#+v[i]/O[1]  
  } else{
    x0 = 0
  }
  
  xdiff = 1
  tol = 1e-6
  iter = 0
  while(xdiff > tol & iter < 1e3){
    x1 = x0-(f(x0, O, tk, Sk)+sk-s[i])/dfdx(x0,O,tk,Sk)
    xdiff = x1-x0
    x0 = x1
    iter = iter+1
    # cat("\niter",iter,"x0:",x0)
  }
  tt[i] = x1
}
res = data.frame(ye = tt)

system.time({
  res$num_fib = fib(mle, pp, s, dt = .1)  
})
system.time({
  res$exp_fib = exp_fib(mle, pp, s, dt=.1)
})


n  = 600
v  = rexp(n)
s  = cumsum(v)
system.time({
  num_fib = fib(mle, pp, s, dt = .001)  
  num_fib = num_fib[num_fib<Tmax]
})
system.time({
  exp_fib = exp_fib(mle, pp, s, dt=.001)
  exp_fib = exp_fib[exp_fib<Tmax]
})
sum((num_fib-exp_fib)^2)

plot(num_fib, exp_fib);abline(0,1)
plot(num_fib-exp_fib)
range(s)
range(integrated_intensity(seq(0,Tmax,1), pp, mle))

plot(res)
sum((res$ye-res$exp_fib)^2)
sum((res$ye-res$num_fib)^2)
