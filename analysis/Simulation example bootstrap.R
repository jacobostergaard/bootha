library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts5.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pars = dowjones$pars
Tmax = dowjones$Tmax
pp   = simulateHawkes(Tmax, pars)

mle  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)


cat(paste0("\nBegin at: ", Sys.time()),"\n")
fn = paste0(datalib,"sim_boot.Rda")
B = 299
ncores = 8
dji.boot = list() 
system.time({
  dji.boot[[1]] = bootstrap_hawkes(pp, mle, Tmax, B, dt = 1, type='fib', parametric=TRUE, ncores=ncores)
})
system.time({
  dji.boot[[2]] = bootstrap_hawkes(pp, mle, Tmax, B, dt = 1, type='fib', parametric=FALSE, ncores=ncores)
})
system.time({
  dji.boot[[3]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=TRUE, ncores=ncores)
})
system.time({
  dji.boot[[4]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=FALSE, ncores=ncores)
})
names(dji.boot) = c("fib_par","fib_npar", "rib_par","rib_npar")

cat(paste0("\nDone at: ", Sys.time()))

# msg = paste0("DJI bootstrap done")
# notify(msg)
# cat("\r",msg, sep="")

par(mfrow=c(2,2))
for(i in 1:4){
  plot(0,0, type='n', xlim=c(0,6500), ylim=c(0,600))
  invisible(lapply(dji.boot[[i]], function(x) lines(x, 1:length(x), type='s', col=add.alpha('black',.3))))
}

save(dji.boot, file=fn)





# MLE  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=TRUE)


x = seq(0,Tmax,1)
L = integrated_intensity(x, pp, mle)
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
mtext("ACF for Simulation Example")
# par(mfrow=c(2,2))

plot_acf(w, lag=100); mtext("Uncorrected waiting times", 3, line=0)
plot_acf(v, lag=100); mtext("Mean corrected waiting times", 3, line=0)
plot_acf(w^2, lag=100); mtext("Uncorrected waiting times squared", 3, line=0)
plot_acf(v^2, lag=100); mtext("Mean corrected waiting times squared", 3, line=0)



# Load bootstrap samples fib
fn = paste0(datalib,"sim_boot.Rda")
load(fn)

lapply(dji.boot,length)

boot_mle <- function(z){
  out = matrix(unlist(lapply(z, function(x) mle(x, Tmax, mle,dt = 0))), nc=3, byrow=TRUE)
  return(out)
}
boot_mle = lapply(dji.boot, boot_mle)

lbls = c("Par FIB", "Nonpar FIB", "Par RIB", "Nonpar RIB")
yLim = matrix(c(0,.03,.5,1,0.01,.05), byrow=TRUE, nc=2)
par(mfrow=c(3,4))
for(i in 1:3){
  for(j in 1:4){
    tmp = boot_mle[[j]][,i]  
    boxplot(tmp, ylim = yLim[i,])
    abline(h=c(mle[i],pars[i]), lty=c(1,3))
    legend("topright", c("mle","sim"), lty=c(1,3), bty='n')
    if(i==1) mtext(lbls[j], side=3, line=2)
  }
}
res   = mle
ses   = sqrt(diag(-solve(exp_hessian(pp, Tmax, mle))))
mleCI = data.frame(sim = pars, mle = res, ci.lo = res-1.96*ses, ci.hi = res+1.96*ses)

pars
boot_ci <- function(x, include_median = FALSE){
  if(include_median){
    bootCI = data.frame(sim = pars, mle = res, t(apply(x, 2, function(y) quantile(y, probs = c(.025,.975,.5)))))
    colnames(bootCI) = c("sim", "mle", "ci.lo", "ci.hi", "median")  
  } else{
    bootCI = data.frame(sim = pars, mle = res, t(apply(x, 2, function(y) quantile(y, probs = c(.025,.975)))))
    colnames(bootCI) = c("sim", "mle", "ci.lo", "ci.hi")  
  }
  return(bootCI)
}
# pars
# round(mle,3)

round(mleCI,4)
lapply(boot_mle, function(x) round(boot_ci(x),4))

set.seed(1234)
x = seq(0,Tmax,.1)
L = integrated_intensity(x, pp, mle)
L = cbind(x,L)

length(pp)
boot1fib <- function(){
  v = rexp(1.5*length(pp))
  out = fib(L,v)
  out = out[out<Tmax]
  return(out)
}
B = 399
system.time({
  bootBfib = replicate(B, boot1fib())  
})

mle_fib = matrix(unlist(lapply(bootBfib,function(x) mle(x, Tmax, mle,dt = 0))),nc=3, byrow=TRUE)

pars
mle
lapply(boot_mle, function(x) apply(x,2, median))
lapply(boot_mle, function(x) apply(x,2, mean))
apply(mle_fib,2, mean)
apply(mle_fib,2, median)
