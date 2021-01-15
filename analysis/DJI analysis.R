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

MLE = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=TRUE)

x = seq(0,Tmax,1)
L = integrated_intensity(x, pp, MLE$mle)
s = interpol(cbind(x,L), x=pp)
w = diff(s)
v = w/mean(w);mean(w);sd(w)

layout(1)

# Plot s vs Poisson process bounds
    # a = qnorm(.975)
    # mu.tmp = length(s)/max(s)
    # mu.hi  = mu.tmp + a*sqrt(mu.tmp)/length(s)
    # mu.lo  = mu.tmp - a*sqrt(mu.tmp)/length(s)
    # plot(s,1:length(s), type='s')
    # abline(0,mu.tmp)
    # abline((mu.hi-mu.tmp)*max(s),mu.tmp, lty=3)
    # abline((mu.lo-mu.tmp)*max(s),mu.tmp, lty=3)


# Plot fitted model
    hi = Tmax+1
    lo = 0
    dt = 365.25*3/12 
    x  = seq(0,Tmax,dt)
    y  = intensity(x, pp, pars)
    z  = pp[pp > lo & pp < hi]
    h  = hist(pp, breaks=seq(0,max(z)+dt,dt), plot=FALSE)
    
    plot(h$mids+dt/2, h$counts, type='l', xlim=c(lo,hi), bty='n', col=add.alpha('steelblue',.5), lwd=5)
    lines(x, y*dt, col=add.alpha('red',.75), lwd=2)
    points(z, rep(0,length(z)), pch=124)


# Plot autocorrelations
    layout(matrix(c(1,1,2,3,4,5), nr=3, nc=2, byrow=TRUE))
    par(mar=c(3,3,1,1), oma=c(0,0,3,0))
    plot_acf(diff(pp), lag=100)
    mtext("ACF for Dow Jones")
    # par(mfrow=c(2,2))
    
    plot_acf(w, lag=100); mtext("Uncorrected waiting times", 3, line=0)
    plot_acf(v, lag=100); mtext("Mean corrected waiting times", 3, line=0)
    plot_acf(w^2, lag=100); mtext("Uncorrected waiting times squared", 3, line=0)
    plot_acf(v^2, lag=100); mtext("Mean corrected waiting times squared", 3, line=0)

    
    
# Load bootstrap samples fib
fn = paste0(datalib,"dji_boot.Rda")
load(fn)

lapply(dji.boot,length)

boot_mle <- function(z){
  out = matrix(unlist(lapply(z, function(x) mle(x, Tmax, MLE$mle,dt = 0))), nc=3, byrow=TRUE)
  return(out)
}
boot_mle = lapply(dji.boot, boot_mle)

yLim = matrix(c(0,.05,.6,1,0,.05), byrow=TRUE, nc=2)
par(mfrow=c(3,4))
for(i in 1:3){
  for(j in 1:4){
      tmp = boot_mle[[j]][,i]  
      boxplot(tmp, ylim = yLim[i,])
      abline(h=c(MLE$mle[i],pars[i]), lty=c(1,3))
  }
}
res   = MLE$mle
ses   = sqrt(diag(-solve(exp_hessian(pp[-(1:194)], Tmax, MLE$mle))))
mleCI = data.frame(mle = res, ci.lo = res-1.96*ses, ci.hi = res+1.96*ses)

pars
boot_ci <- function(x, include_median = FALSE){
  if(include_median){
    bootCI = data.frame(mle = res, t(apply(x, 2, function(y) quantile(y, probs = c(.025,.975,.5)))))
    colnames(bootCI) = c("mle", "ci.lo", "ci.hi", "median")  
  } else{
    bootCI = data.frame(mle = res, t(apply(x, 2, function(y) quantile(y, probs = c(.025,.975)))))
    colnames(bootCI) = c("mle", "ci.lo", "ci.hi")  
  }
  return(bootCI)
}

round(mleCI,4)
lapply(boot_mle, function(x) round(boot_ci(x),4))


# booter <- function(B, type, parametric, dt=1, ncores=6){
#   tmp = bootstrap_hawkes(pp, pars, Tmax, B, type, parametric,  dt, ncores)
#   cat("    Bootstrap done, estimating mle...")
#   res = parallel::mclapply(tmp, function(pp) mle(pp, Tmax, pars, dt=0), mc.cores=ncores)
#   # res = lapply(tmp, function(pp) mle(pp, Tmax, pars) )
#   res = matrix(unlist(res),nc=3, byrow=TRUE)
#   colnames(res) = c("mu","n","B")
#   return(res)
# }
    
# dt = 1  => 399 rib samples in about  5h
# dt = 5  => 399 rib samples in about  1h
# dt = 10 => 399 rib samples in about .5h
    
