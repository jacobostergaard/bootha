# Analysis of DJI
library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"

savePDF = FALSE
savePDF = TRUE

# Make point process of Dow Jones Index log return exceedances
fn  = paste0(datalib,"DJI.csv")
dji = read.csv(fn)
dji$Date = as.Date(dji$Date)
dji$log.return = c(NA,diff(log(dji$Close)))
dji = dji[,c(1,5,7,8)] # relevant columns

date.from  = as.Date("1994-01-01")
date.to    = as.Date("2010-12-31")
idx        = dji$Date>date.from & dji$Date < date.to # same period as in Embrechts paper
names(dji) = c("date", "close", "volume", "log.return")
dji        = dji[idx,]

qts = quantile(dji$log.return, probs = c(.1,.9), na.rm = TRUE) 
idx = which(dji$log.return < qts[1] | dji$log.return > qts[2]) # find days where the log return exceeds the .1 and .9 quantiles

x    = as.numeric(dji$date[idx])
x    = x-min(x)+1
y    = dji$log.return[idx]

pars = pars0 = c(0.018,0.74,0.021) # Parameters from the Embrechts paper
Tmax0 = max(x)
dji  = list(pp.neg = as.numeric(x[y<0]), pp.pos = as.numeric(x[y>0]))

# Estimate MLE and rescaled time transformed waiting times
pp0   = dji$pp.neg # choose events of negative exceedances

t_norm = 1
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle1  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

t_norm = 5
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle5  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

t_norm = 10
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle10  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

t_norm = 15
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle15  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

t_norm = 20
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle20  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

t_norm = 25
# t_norm = max(x)/length(pp)
pp = pp0/t_norm
Tmax = Tmax0/t_norm
pars[c(1,3)]=pars0[c(1,3)]*t_norm
mle25  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

Tmax0/Tmax

mle5/mle1
mle10/mle1
mle15/mle1
mle20/mle1
mle25/mle1
Tmax0/Tmax

Tmax*mle25
Tmax0*mle1
