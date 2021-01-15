library(emhawkes)
library(hawkes)
library(ppboot)
library(misc)
library(parallel)
clean_up()

tic = Sys.time()

B    = 100
mu   = 1
a    = 4
b    = 5
Tmax = 50

spec <- new("hspec", mu=mu, alpha=a, beta=b)

set.seed(1234)
obs = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
fit = quiet(emhawkes::hfit(spec, diff(c(0,obs))))

theta0 = c(mu,a,b)
thetah = fit$estimate
llhobs = fit$maximum

mle = c(thetah,llhobs)
names(mle) = c("mu","alpha","beta","loglik")

boot <- function(obs, B, pars, parametric = TRUE, type='fib', dt=1e-2){
  
  mu.in = pars[1]
  a.in = pars[2]
  b.in = pars[3]
  
  x = seq(0,Tmax,dt)
  l = data.frame(x=x,y=l_exp_kernel(x,obs,mu.in,a.in,b.in))
  L = data.frame(x=x,y=L_exp_kernel(x,obs,mu.in,a.in,b.in))
  s = interpol(as.matrix(L),x=obs)
  if(type=='fib'){
    kern = NULL
  } else{
    kern <- function(x,y){
      if(is.null(y)){
        y = numeric(0)
      }
      out = L_exp_kernel(x,y,mu.in,a.in,b.in)
    }  
  }
  if(parametric){
    s.events = NULL
  } else{
    s.events = s
  }
  
  Smax = max(L$y)
  bb = pp.boot(L = L, Smax = Smax, kernel = kern, B = B, s.events = s.events, parallel = TRUE)
 
  return(bb) 
}

fit_boot <- function(bb, pars){
  spec = new("hspec", mu=pars[1], alpha=pars[2], beta=pars[3])
  tmp.fit = lapply(bb, function(x) fit = quiet(emhawkes::hfit(spec, diff(c(0,x)), lambda0=pars[1])))
  tmp.est = matrix(unlist(lapply(tmp.fit, function(x) x$estimate)), nc=3, byrow=TRUE)
  tmp.est = cbind(tmp.est,unlist(lapply(tmp.fit, function(x) x$maximum)))
  colnames(tmp.est) = c("mu","alpha","beta","loglik")
  rownames(tmp.est) = 1:length(bb)
  return(tmp.est)  
}

# FIB
  cat("\nRestricted parametric FIB")
  system.time({ parfib.r = boot(obs,B=B, theta0, parametric = TRUE, type='fib') })
  cat("\nUnrestricted parametric FIB")
  system.time({ parfib.u = boot(obs,B=B, thetah, parametric = TRUE, type='fib') })
  cat("\nRestricted nonparametric FIB")
  system.time({ nparfib.r = boot(obs,B=B, theta0, parametric = FALSE, type='fib') })
  cat("\nUnrestricted nonparametric FIB")
  system.time({ nparfib.u = boot(obs,B=B, thetah, parametric = FALSE, type='fib') })
    
# RIB
  cat("\nRestricted parametric RIB")
  system.time({ parrib.r = boot(obs,B=B, theta0, parametric = TRUE, type='rib') })
  cat("\nUnrestricted parametric RIB")
  system.time({ parrib.u = boot(obs,B=B, thetah, parametric = TRUE, type='rib') })
  cat("\nRestricted nonparametric RIB")
  system.time({ nparrib.r = boot(obs,B=B, theta0, parametric = FALSE, type='rib') })
  cat("\nUnrestricted nonparametric RIB")
  system.time({ nparrib.u = boot(obs,B=B, thetah, parametric = FALSE, type='rib') })

  print("MLE:")
  print(round(mle,3))

  cat("\nFitting: Restricted parametric FIB")
  system.time({ (fit.parfib.r = fit_boot(parfib.r, theta0)) })
  cat("\nFitting: Unrestricted parametric FIB")
  system.time({ (fit.parfib.u = fit_boot(parfib.u, thetah)) })
  cat("\nFitting: Restricted nonparametric FIB")
  system.time({ (fit.nparfib.r = fit_boot(nparfib.r, theta0)) })
  cat("\nFitting: Unrestricted nonparametric FIB")
  system.time({ (fit.nparfib.u = fit_boot(nparfib.u, thetah)) })
    
  cat("\nFitting: Restricted parametric RIB")
  system.time({ (fit.parrib.r = fit_boot(parrib.r, theta0)) })
  cat("\nFitting: Unrestricted parametric RIB")
  system.time({ (fit.parrib.u = fit_boot(parrib.u, thetah)) })
  cat("\nFitting: Restricted nonparametric RIB")
  system.time({ (fit.nparrib.r = fit_boot(nparrib.r, theta0)) })
  cat("\nFitting: Unrestricted nonparametric RIB")
  system.time({ (fit.nparrib.u = fit_boot(nparrib.u, thetah)) })


par(mfrow=c(2,1), mar=c(3,8,1,1), oma=c(1,1,1,1))
col1 = add.alpha('red',.75)
col2 = add.alpha('dodgerblue',.75)
col3 = add.alpha('orange',.75)
col4 = add.alpha('green',.75)
plot(obs, 1:length(obs), type='l', lwd=2, ylim=c(0,4*mu/(1-a/b)*Tmax))
for(i in 1:B){
  tmp.obs = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  lines(tmp.obs,1:length(tmp.obs), col=add.alpha('black',.35), lty=1, lwd=2)
  
  tmp.obs = parfib.r[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col1, lty=1)
  
  tmp.obs = nparfib.r[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col1, lty=2)
  
  tmp.obs = parfib.u[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col2, lty=1)
  
  tmp.obs = nparfib.u[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col2, lty=2)
  
  tmp.obs = parrib.r[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col3, lty=1)
  
  tmp.obs = nparrib.r[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col3, lty=2)
  
  tmp.obs = parrib.u[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col4, lty=1)
  
  tmp.obs = nparrib.u[[i]]
  lines(tmp.obs,1:length(tmp.obs), col=col4, lty=2)
}
legend("topleft",c("Obs","Hawkes","Restricted FIB","Un-restricted FIB","Restricted RIB","Un-restricted RIB","Parametric","Non-parametric"), col=c('black',add.alpha('black',.5),col1,col2,col3,col4,1,1), lty=c(1,1,1,1,1,1,1,2), lwd=c(2,2,1,1,1,1,1,1), bty='n')

cexs=2
plot(obs, rep(0,length(obs)), pch='.', ylim=c(0,9), yaxt='n', ylab="", cex=cexs)
for(i in 1:B){
  tmp.obs = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  points(tmp.obs,rep(0+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = parfib.r[[i]]
  points(tmp.obs,rep(1+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = nparfib.r[[i]]
  points(tmp.obs,rep(2+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = parfib.u[[i]]
  points(tmp.obs,rep(3+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = nparfib.u[[i]]
  points(tmp.obs,rep(4+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = parrib.r[[i]]
  points(tmp.obs,rep(5+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = nparrib.r[[i]]
  points(tmp.obs,rep(6+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = parrib.u[[i]]
  points(tmp.obs,rep(7+i/B,length(tmp.obs)), pch='.', cex=cexs)
  
  tmp.obs = nparrib.u[[i]]
  points(tmp.obs,rep(8+i/B,length(tmp.obs)), pch='.', cex=cexs)
}

lbls = c("Observed", "Simulated Hawkes", "Res. par. FIB", "Res. nonpar. FIB", "Unres. par. FIB", "Unres. nonpar. FIB", "Res. par. RIB", "Res. nonpar. RIB", "Unres. par. RIB", "Unres. nonpar. RIB")
tcks = c(0,0:8+2/3)
axis(2, tcks, lbls, las=1, cex=.3)
abline(h=1/6+0:8, lty=3)

toc = Sys.time()
cat("\nTotal time:",as.numeric(difftime(toc,tic, units = "secs")))
