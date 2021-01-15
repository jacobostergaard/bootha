library(emhawkes)
library(hawkes)
library(ppboot)
library(misc)
clean_up()


mu    = 1
a    = 4
b    = 5
Tmax = 100

spec <- new("hspec", mu=mu, alpha=a, beta=b)

set.seed(1234)
hobs = c(0,unlist(hawkes::simulateHawkes(mu,a,b,Tmax)))

hawkes::likelihoodHawkes(mu,a,b,hobs)

res = quiet(emhawkes::hfit(spec, diff(hobs)))

theta.tru = c(mu,a,b)
theta.hat = res$estimate

theta = theta.tru
H = hobs

LRT <- function(theta.1, theta.0, H){
  
  theta = theta.0
  l0 = -hawkes::likelihoodHawkes(theta[1],theta[2],theta[3],H)  
  theta = theta.1
  l1 = -hawkes::likelihoodHawkes(theta[1],theta[2],theta[3],H)  
  
  out = -2*(l0-l1)
  
  return(out)
}

LRT.obs = LRT(theta.hat, theta.tru, hobs)


mu.in = mu#mle$estimate[1]
a.in = a#mle$estimate[2]
b.in = b#mle$estimate[3]
bootspec <- new("hspec", mu=mu.in, alpha=a.in, beta=b.in)

dt = 1e-2
x = seq(0,Tmax,dt)
l = data.frame(x=x,y=l_exp_kernel(x,hobs,mu.in,a.in,b.in))
L = data.frame(x=x,y=L_exp_kernel(x,hobs,mu.in,a.in,b.in))

s  = interpol(as.matrix(L),x=hobs)

kern <- function(x,y){
  if(is.null(y)){
    y = numeric(0)
  }
  out = L_exp_kernel(x,y,mu.in,a.in,b.in)
}

B    = 399
Smax = max(L$y)

library(parallel)
system.time({
  bb = pp.boot(L,Smax,kern,B,NULL,parallel = TRUE)
  # bb = pp.boot(L,Smax,NULL,B,NULL,parallel = TRUE)
})


f <- function(i){
  bb.tmp = c(0,bb[[i]])
  res.tmp = quiet(emhawkes::hfit(spec, diff(bb.tmp)))
  lrt1 = LRT(res.tmp$estimate, theta.tru, bb.tmp)
  lrt2 = LRT(res.tmp$estimate, theta.hat, bb.tmp)
  out = c(lrt1, lrt2, res.tmp$estimate)
  names(out) = c("lrt.0","lrt.hat","mu","alpha","beta")
  return(out)
}

system.time({
  Bres = parallel::mclapply(1:B,f,mc.cores = 4)  
})

LRT.boot1 = matrix(unlist(Bres),nc=5, byrow = TRUE)[,1]
LRT.boot2 = matrix(unlist(Bres),nc=5, byrow = TRUE)[,2]

col1 = add.alpha('dodgerblue',.5)
col2 = add.alpha('red',.5)
col3 = add.alpha('black',.5)
hist(LRT.boot1, breaks=seq(0,50,1), border=NA, col=col1, prob=TRUE, ylim=c(0,.2), xlim=c(0,100))
par(new=TRUE)
hist(LRT.boot2, breaks=seq(0,500,1), border=NA, col=col2, prob=TRUE, ylim=c(0,.2), xlim=c(0,100))
abline(v=LRT.obs)
curve(dchisq(x,3), add=TRUE, col=col3, lwd=2, lty=3)

