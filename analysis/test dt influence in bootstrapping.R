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

pp.embrechts = pp

set.seed(12)
Tmax
pp.simulated = simulateHawkes(Tmax, pars)

length(pp.simulated) # simulated observations
length(pp.embrechts)

(mle.obs = mle(pp.embrechts, Tmax, pars))
log_lik(pp.embrechts, Tmax, pars)
log_lik(pp.embrechts, Tmax, mle.obs)

(mle.sim = mle(pp.simulated, Tmax, pars))
log_lik(pp.simulated, Tmax, pars)
log_lik(pp.simulated, Tmax, mle.sim)


booter <- function(B, type, parametric, dt=1, ncores=8){
  tmp =  bootstrap_hawkes(pp, pars, Tmax, B, type, parametric,  dt, ncores)
  res = lapply(tmp, function(pp) mle(pp, Tmax, pars))
  res = matrix(unlist(res),nc=3, byrow=TRUE)
  colnames(res) = c("mu","n","B")
  return(res)
}

B = 5
# system.time({
#   bootstrap_hawkes(pp.embrechts, pars, Tmax, B, type='fib', parametric=TRUE,  dt=10, ncores=5)
# })
# system.time({
#   bootstrap_hawkes(pp.simulated, pars, Tmax, B, type='fib', parametric=TRUE,  dt=10, ncores=5)
# })
# 
# system.time({
#   bootstrap_hawkes(pp.embrechts, pars, Tmax, B, type='fib', parametric=TRUE,  dt=10, ncores=1)
# })
# system.time({
#   bootstrap_hawkes(pp.simulated, pars, Tmax, B, type='fib', parametric=TRUE,  dt=10, ncores=1)
# })

B = 299
T.tmp =3000

N1 = pp.embrechts[pp.embrechts<T.tmp]
N2 = pp.simulated[pp.simulated<T.tmp]
length(N1)
length(N2)


system.time({
  res1.N1 = bootstrap_hawkes(N1, mle.obs, T.tmp, B, type='rib', parametric=TRUE,  dt=1, ncores=6)
})
system.time({
  res5.N1 = bootstrap_hawkes(N1, mle.obs, T.tmp, B, type='rib', parametric=TRUE,  dt=5, ncores=6)
})

system.time({
  res1.N2 = bootstrap_hawkes(N2, mle.sim, T.tmp, B, type='rib', parametric=TRUE,  dt=1, ncores=6)
})
system.time({
  res5.N2 =bootstrap_hawkes(N2, mle.sim, T.tmp, B, type='rib', parametric=TRUE,  dt=5, ncores=6)
})


plot(N1, 1:length(N1), type='s')
lines(N2, 1:length(N2), type='s', col='red')


make_mle <- function(tmp){
  res = parallel::mclapply(tmp, function(pp) mle(pp, T.tmp, pars), mc.cores=6)
  res = matrix(unlist(res),nc=3, byrow=TRUE)
  colnames(res) = c("mu","n","B")
  return(res)  
}

fn = paste0(datalib,"dt_test_boot.Rda")
dtboot = list("dt1.N1"=res1.N1, "dt5.N1"=res5.N1, "dt1.N2"=res1.N2, "dt5.N2"=res5.N2)
save(dtboot, file=fn)

mle5.N1 = parallel::mclapply(res5.N1, function(pp) mle(pp, T.tmp, pars), mc.cores = 6)

mle1.N1 = make_mle(res1.N1)
mle1.N2 = make_mle(res1.N2)
mle5.N1 = make_mle(res5.N1)
mle5.N2 = make_mle(res5.N2)

fn = paste0(datalib,"dt_test_mle.Rda")
boot = list("dt1.N1"=mle1.N1, "dt5.N1"=mle5.N1, "dt1.N2"=mle1.N2, "dt5.N2"=mle5.N2)
save(boot, file=fn)

parnum=1
parname=expression(mu)
cols = add.alpha(c('steelblue', 'tomato2', 'tomato2'),.8)
cols = add.alpha(c('steelblue','tomato2'),.5)
plotpar <- function(parnum, parname="", cols=add.alpha('black',.5), yl=NULL){
  tmp = numeric(0)
  for(i in 1:4){
    tmp = cbind(tmp,boot[[i]][,parnum])  
  }
  if(is.null(yl)){
    yl = c(0,max(tmp))  
  }
  
  boxplot(tmp, col=cols[1], border='black', pch=16, cex=.75, xaxt='n', ylim=yl)
  abline(h=pars[parnum], lty=3, lwd=2)
  # abline(h=res[parnum], lty=1, lwd=2)
  tmp = apply(tmp,2, make_ci)
  for(i in 1:4){
    # segments(i,tmp[1,i],i,tmp[2,i], lwd=50, col=cols[2], lend='butt')
  }
  lbls = c("Embrechts dt=1", "Embrechts dt=5", "Simulated dt=1", "Simulated dt=5")
  text(1:4,.45*par("usr")[3], lbls, srt=0, xpd = NA, cex=.85)
  mtext(parname,side = 3, line=0)
}

layout(1:3)
par(mar=c(3,2,2,1), bty='n', oma=c(1,1,6,1))
plotpar(1, expression(mu), cols, yl=c(0,.05))
plotpar(2, "n", cols, yl=c(0,1))
plotpar(3, expression(beta), cols, yl=c(0,.05))
lgnd = c("MLE", "Sim. params", "Int-Q range","95% CI")
legend("topleft",lgnd, col=c('black','black',cols), lty=c(1,3,NA,NA), pch=c(NA,NA,15,15), cex=1, bty='n', pt.cex = 2)
mtext("Bootstrap Embrechts", 3, line=3, outer=TRUE, cex=1.5)  
mtext("Using 299 bootstrap samples", 3, line=1.5, outer=TRUE)  
