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


make_cmp <- function(x,cut){
  xmax  = 0
  addon = 0
  cut   = c(cut, max(x))
  Z.cmp = 0

  for(i in 1:length(cut)){
    idx   = which( x >= xmax & x <= cut[i])
    x.tmp = x[idx]
    Z.tmp = integrated_intensity(x.tmp,pp,pars)
    Z.tmp = Z.tmp+addon
    Z.cmp = c(Z.cmp,Z.tmp[-1])
    addon = max(Z.tmp)
    xmax  = max(x.tmp)  
  }
  return(Z.cmp)
}


dt  = .01
x0  = seq(0,100,dt)
Z0 = integrated_intensity(x0,pp,pars) # finegrained compensator
x   = sort(c(0,sample(x0[c(-1,-length(x))], size = 300, replace = FALSE),100)) # sample compensator at different points
cut = sort(c(50.1,80,99,24))
Z = integrated_intensity(x,pp,pars) # sampled compensator

layout(1)
plot(x0,Z0, type='l') # plot "true" compensator
points(x,Z, col='red', pch=16, cex=.5) # plot sampled compensator

sum((Z0-make_cmp(x0,cut))^2) # MSE with "true" compensator
sum((Z-make_cmp(x,cut))^2)   # MSE with sampled compensator


