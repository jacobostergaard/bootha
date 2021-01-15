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

f    <- function(x,O,Sm,v) O[1]*x+O[2]*(1-exp(-O[3]*x))*Sm-v
dfdx <- function(x,O,Sm) O[1]+O[2]*O[3]*exp(-O[3]*x)*Sm 
ye   <- function(O, v){
  
  n  = length(v)
  w  = numeric(n)
  Sm = 1
  
  i=1
  w[i] = v[i]/O[1]
  for(i in 2:n){
    x0 = 0
    xdiff = 1
    tol = 1e-6
    iter = 0
    while(xdiff > tol & iter < 10){
      x1 = x0-f(x0,O,Sm,v[i])/dfdx(x0,O,Sm)
      xdiff = x1-x0
      x0 = x1
      iter = iter+1
      # cat("\niter",iter,"x0:",x0)
    }
    # w[i] = uniroot(f, c(0,100))$root
    w[i] = x0
    Sm = exp(-O[3]*w[i])*Sm+1
  }
  
  out = cumsum(w)
  return(out)
}


# g <- function(x) f(x)^2
# nlm(g, 0)

# x = seq(0,10,.1)
# plot(x, f(x), type='l')


n = 100
v  = rexp(n)

system.time({
  num_rib = rib(mle, cumsum(v), dt=1)  
})
system.time({
  ye_rib = ye(mle,v)
})
system.time({
  exp_rib = exp_rib(mle, v)
})

res = data.frame(ye=ye_rib,num=num_rib, exp = exp_rib)


sum((res$ye-res$num)^2)
sum((res$ye-res$exp)^2)
sum((res$exp-res$num)^2)


par(mar=c(4,4,1,1))
plot(res)
abline(0,1)
