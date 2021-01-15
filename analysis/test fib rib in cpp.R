library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

dji <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/DJI.csv")
dji$Date = as.Date(dji$Date)
dji$log.return = c(NA,diff(log(dji$Close)))

date.from = as.Date("1994-01-01")
date.to = as.Date("2010-12-31")

idx.dji = dji$Date>date.from & dji$Date < date.to
idx = c(1,5,7,8)
dji.tmp = dji[idx.dji,idx]
pfc = dji.tmp
names(pfc) = c("date", "close.dji", "volume.dji", "log.return.dji")

qts = quantile(pfc$log.return, probs = c(.01,.99), na.rm = TRUE)
qts = quantile(pfc$log.return, probs = c(.1,.9), na.rm = TRUE)

idx = which(pfc$log.return < qts[1] | pfc$log.return > qts[2])

x = as.numeric(pfc$date[idx])
x = x-min(x)+1
y = pfc$log.return[idx]

mu = 0.018
n = .74
d = 0.021
pp = as.numeric(x[y<0])

pars = c(mu,n,d)

dt = .5
x   = seq(0,6144,dt)
z = intensity(x,pp,pars)
Z = integrated_intensity(x,pp,pars)

# plot(x,z, type='l')
# plot(x,dt*cumsum(z), type='l')
# lines(x,Z, col='red')

Tmax = max(pp)

hawkes::likelihoodHawkes(mu, n*d, d, pp)
criticality::log_lik(pp,Tmax, mu, n, d)
log_lik(pp, Tmax, pars, 0)
log_lik(pp, Tmax, pars)
log_lik(pp, Tmax, pars, 1)

ppy = interpol(cbind(x,Z), x=pp)

s  = cumsum(rexp(length(pp))) # sample events in S-time

system.time({
  b.fib = fib_cpp(pars, pp, s, dt=.5)
  # b.fib = fib_cpp(cbind(x,Z),s)
})
system.time({
   b.rib = rib_cpp(pars,s, dt=.5)  
})

Z.rib = integrated_intensity(x, b.rib, pars)

col1 = add.alpha('black',.95)
col2 = add.alpha('dodgerblue',.95)
col3 = add.alpha('tomato2',.75)

# plot(x,Z, col=col1, type='l', xlim=c(0,1000), ylim=c(0,35))
# lines(x,Z.rib, col=col3)
# segments(pp, 0, pp, ppy, lty=1, col=col1); segments(0, ppy, pp, ppy, lty=1, col=col1)
# segments(b.fib, 0, b.fib, s, lty=1, col=col2); segments(0, s, b.fib, s, lty=1, col=col2)
# segments(b.rib, 0, b.rib, s, lty=1, col=col3); segments(0, s, b.rib, s, lty=1, col=col3)

par(mar=c(1,3,0,0), oma=c(3,0,0,0), bty='n')
layout(1:2, height = c(5,1))
plot(x,Z, col=col1, type='l', xaxt='n')#, xlim=c(0,1000), ylim=c(0,35))
lines(x,Z.rib, col=col3)
plot(0,0, type='n', xlim=c(0,6000), ylim=c(0,4), yaxt='n')
points(pp, rep(1,length(pp)), pch=124, col=col1)
points(b.fib, rep(2,length(b.fib)), pch=124, col=col2)
points(b.rib, rep(3,length(b.rib)), pch=124, col=col3)
# points(rep(0,length(ppy)), ppy, pch=95, col=col1)
# points(rep(0,length(s)), s, pch=95, col=col2)
length(b.rib)
length(b.fib)
length(pp)

f = function(pars){
  -log_lik(pp, Tmax, pars)
}
system.time({
  res = constrOptim(pars, f, ui = diag(3), ci=rep(0,3), method="Nelder-Mead")
})
res$par
pars

g = function(pars){
  -log_lik(b.rib, Tmax, pars)
}

system.time({
  res = constrOptim(pars, g, ui = diag(3), ci=rep(0,3), method="Nelder-Mead")
})
res$par
pars


# system.time({
#   b.rib = rib_cpp(pars,s[1:50], dt=1)  
# })
xt = c(50, 100, 200, 400)
t1 = c(0.05, 0.37, 1.88, 12.78) # dt=.5
t2 = c(0.03, 0.17, 1.07, 6.82)  # dt=1

# scales approx linear with dt, i.e. time*dt
# scales arppx quadratic in #s, i.e. (#s)^2
# implies: approx. time = (#s)^2*dt

layout(1)
plot(xt, sqrt(t1), col=col2, type='l')
lines(xt, sqrt(t2), col=col3)

xt  = c(50, 100, 200, 400)
dts = c(.25,.5,1)
rtime = matrix(nr=length(xt), nc=length(dts))
for(i in 1:length(xt)){
  for(j in 1:length(dts)){
    cat(i,j, "\n")
    rtime[i,j] = system.time({
                                b.tmp = rib_cpp(pars,s[1:xt[i]], dt=dts[j])
                              })[3]
  }
}

colnames(rtime) = paste0("dt",dts)
rownames(rtime) = paste0("t",xt)
rtime

layout(1)
plot(0,0, type='n', xlim=c(0,max(xt)), ylim=sqrt(c(0,max(rtime))), yaxt='n')
for(i in 1:length(dts)){
  lines(xt,sqrt(rtime[,i]), col=i)
}
legend("topleft",legend = dts, fill=1:length(dts), bty='n')
yx = c(seq(0,1,.1),2:10)
axis(2, at=sqrt(yx), labels=paste0(yx,"s"), las=1, cex.axis=.5)



plot(pp,b.fib, pch=16, cex=.5)
points(pp,b.rib, pch=16, cex=.5, col='red')
points(b.fib,b.rib, pch=16, cex=.5, col='blue')




g = function(pars){
  # -log_lik(b.rib, Tmax, pars)
  # -log_lik(b.fib, Tmax, pars)
  -log_lik(pp, Tmax, pars)
}

n0  = .9
mu0 = n0*length(b.rib)/Tmax
B0  = .021
ui0 = diag(c(1,-1,1))
system.time({
  res = constrOptim(c(mu0, .5, .01), g, ui = ui0, ci=c(0,-1,0), method="Nelder-Mead")
})
res$par
res$par[2]*res$par[3]
pars

-log_lik(b.rib, Tmax, pars)
-log_lik(b.fib, Tmax, pars)
-log_lik(pp, Tmax, pars)

