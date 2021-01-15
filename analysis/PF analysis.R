library(ppboot)
library(misc)
misc::clean_up()

dji <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/DJI.csv")
ndx <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/NDX.csv")
spc <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/SPC.csv")


dji$Date = as.Date(dji$Date)
ndx$Date = as.Date(ndx$Date)
spc$Date = as.Date(spc$Date)

range(dji$Date)
range(ndx$Date)
range(spc$Date)

dji.col = add.alpha('steelblue',.75)
ndx.col = add.alpha('tomato2',.75)
spc.col = add.alpha('darkgoldenrod2',.75)

plot(spc$Date, spc$Close, type='l', col=spc.col, lwd=2, ylim=c(0,3e4))
lines(ndx$Date, ndx$Close, col=ndx.col, lwd=2)
lines(dji$Date, dji$Close, col=dji.col, lwd=2)

spc$log.return = c(NA,diff(log(spc$Close)))
ndx$log.return = c(NA,diff(log(ndx$Close)))
dji$log.return = c(NA,diff(log(dji$Close)))

plot(spc$Date, spc$log.return, type='l', col=spc.col, lwd=2)
lines(ndx$Date, ndx$log.return, col=ndx.col, lwd=2)
lines(dji$Date, dji$log.return, col=dji.col, lwd=2)


date.from = as.Date("1997-10-01")
date.to = as.Date("2010-03-01")

idx.dji = dji$Date>date.from & dji$Date < date.to
idx.ndx = ndx$Date>date.from & ndx$Date < date.to
idx.spc = spc$Date>date.from & spc$Date < date.to
head(dji)
idx = c(1,5,7,8)
dji.tmp = dji[idx.dji,idx]
ndx.tmp = ndx[idx.ndx,idx]
spc.tmp = spc[idx.spc,idx]
pfc = dplyr::left_join(dji.tmp, ndx.tmp, by = "Date")
pfc = dplyr::left_join(pfc, spc.tmp, by = "Date")

names(pfc) = c("date", "close.dji", "volume.dji", "log.return.dji", "close.ndx", "volume.ndx", "log.return.ndx", "close.spc", "volume.spc", "log.return.spc")
pfc$log.return = pfc$log.return.dji+pfc$log.return.ndx+pfc$log.return.spc
head(pfc)


plot(pfc$date, pfc$log.return, type='l')

qts = quantile(pfc$log.return, probs = c(.01,.99), na.rm = TRUE)
qts = quantile(pfc$log.return, probs = c(.1,.9), na.rm = TRUE)
abline(h=qts, lty=3)

idx = which(pfc$log.return < qts[1] | pfc$log.return > qts[2])

x = pfc$date[idx]
y = pfc$log.return[idx]

par(bty='n', mar=c(2,2,0,0), oma=c(0,0,0,0), mfrow=c(2,1))

plot(x[y>0], y[y>0], type='l')
plot(x[y<0], y[y<0], type='l')

z = as.numeric(x)
# z = (z-min(z))/365.25
z = (z-min(z))

plot(z, rep(1, length(z)), pch=124)

# tmp = as.data.frame(table(z))
# tmp = tmp[order(tmp$Freq,decreasing=TRUE),]

layout(1)
min(x)
# h = hist(z, breaks=seq(0,13,3/12), border=NA, col=dji.col)
h = hist(z, breaks=seq(0,max(z+91), 365.25*3/12), border=NA, col=dji.col)

Tmax = ceiling(max(z))
pp = z[z<Tmax]

mu0 = length(pp)/Tmax
n0  = 1-(1/mu0)
n0  = 1-mu0
B0  = 1

mu0 = .008
B0  = .017
n0  = .914

dt = 1
x = seq(0,Tmax,dt)

pars = c(mu0,n0,B0)
f = function(pars){
  -criticality::log_lik(pp, Tmax, pars[1], pars[2], pars[3])
}

system.time({
  res = constrOptim(c(mu0,n0,B0), f, ui = diag(3), ci=rep(0,3), method="Nelder-Mead")
})


pars = res$par

pars
pars = c(0.008, 0.914, 0.017)

L = data.frame(x=x,y=L_exp_kernel(x, pp, pars[1], pars[2], pars[3]))
s = interpol(as.matrix(L), x=pp)

hist(diff(s), breaks=100, prob=TRUE)
curve(dexp, add=TRUE, from=1e-8, col=add.alpha('red',.75), lwd=2)


# Plot fitted model

hi = Tmax+1
lo = 0
dt = 365.25*3/12 
x = seq(0,Tmax,dt)
y = criticality::cond_int(x, pp, pars[1], pars[2], pars[3])
z = pp[pp > lo & pp < hi]

h = hist(pp, breaks=seq(0,max(z)+90,365.25*3/12), plot=FALSE)

set_cobalt_plot_bg()
plot(h$mids, h$counts, type='l', xlim=c(lo,hi), bty='n', col=dji.col, lwd=5)
lines(x, y*dt, col=add.alpha('red',.75), lwd=2)
points(z, rep(0,length(z)), pch=124)


plot(x, L_exp_kernel(x, pp, pars[1], pars[2], pars[3]))




